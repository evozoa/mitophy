"""Stage 2: stream-parse RefSeq mitochondrial GenBank records; write genome stats table and per-marker protein FASTA.

Marker assignment is name-based: CDS /gene and /product qualifiers are normalised and matched against synonym lists in
config/markers.yml. Protein sequences come from /translation (never re-translated) unless it is missing, in which case
the CDS is translated with its /transl_table (and /codon_start). Unmatched product names are logged for curation.
"""
from __future__ import annotations
import csv, gzip, logging, re
from collections import Counter, defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from .config import Config
from .cache import hash_inputs, up_to_date, write_manifest
from .taxonomy import GroupMapper
from . import fetch_mito

log = logging.getLogger("mitophy")

GENOME_COLUMNS = ["accession", "organism", "taxid", "group", "supergroup", "lineage", "definition", "length", "gc",
                  "n_cds", "n_trna", "n_rrna", "transl_table", "n_markers", "markers", "topology", "date"]

_ws = re.compile(r"\s+")
_paren = re.compile(r"\s*\([^)]*\)")

def normalise(s: str) -> str:
    s = s.strip().lower().replace("_", " ").replace("-", " ")
    s = _ws.sub(" ", s)
    s = s.rstrip(".;:, ")
    if s.startswith("mt "):
        s = s[3:]
    return s

def _variants(s: str) -> list[str]:
    n = normalise(s)
    out = [n]
    n2 = normalise(_paren.sub("", s))
    if n2 != n:
        out.append(n2)
    for suf in (" protein", " gene product", " precursor", " like", " homolog"):
        if n.endswith(suf):
            out.append(n[: -len(suf)])
    return out

class MarkerMatcher:
    def __init__(self, markers_cfg: dict):
        self.lookup: dict[str, str] = {}
        for m in markers_cfg["markers"]:
            for syn in [m["id"], *m.get("synonyms", [])]:
                for v in _variants(str(syn)):
                    self.lookup.setdefault(v, m["id"])
        self.ignore = [normalise(x) for x in markers_cfg.get("ignore_products", [])]
        self.ids = [m["id"] for m in markers_cfg["markers"]]

    def match(self, gene: str | None, product: str | None) -> str | None:
        for q in (gene, product):
            if not q:
                continue
            for v in _variants(q):
                if v in self.lookup:
                    return self.lookup[v]
        return None

    def is_ignored(self, product: str | None) -> bool:
        if not product:
            return True
        n = normalise(product)
        return any(x in n for x in self.ignore)


def _taxid(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for x in f.qualifiers.get("db_xref", []):
                if x.startswith("taxon:"):
                    return x.split(":", 1)[1]
    return ""

def _protein(feat, rec) -> str | None:
    tr = feat.qualifiers.get("translation")
    if tr:
        return tr[0].replace("*", "")
    try:
        table = int(feat.qualifiers.get("transl_table", ["1"])[0])
        nt = feat.extract(rec.seq)
        cs = int(feat.qualifiers.get("codon_start", ["1"])[0]) - 1
        nt = nt[cs:]
        nt = nt[: len(nt) - len(nt) % 3]
        aa = str(Seq(nt).translate(table=table)).rstrip("*")
        if "*" in aa or len(aa) < 20:
            return None
        return aa
    except Exception:
        return None


def parse_records(gbff_paths: list[Path]):
    for p in gbff_paths:
        opener = gzip.open if str(p).endswith(".gz") else open
        with opener(p, "rt") as fh:
            yield from SeqIO.parse(fh, "genbank")


def run(cfg: Config) -> Path:
    outdir = cfg.sub("mito")
    gbffs = fetch_mito.gbff_files(cfg)
    ih = hash_inputs(files=gbffs + [cfg.root / "config" / "markers.yml", cfg.root / "data" / "taxonomy" / "supergroups.yml"],
                     objs=[cfg.pipeline.get("filters", {})])
    genomes_tsv = outdir / "genomes.tsv"
    if up_to_date(outdir, ih, [genomes_tsv], cfg.force):
        log.info("extract-mito up to date")
        return genomes_tsv

    matcher = MarkerMatcher(cfg.markers)
    mapper = GroupMapper(cfg.supergroups)
    filt = cfg.pipeline.get("filters", {}).get("mito", {})
    excl = re.compile(filt.get("exclude_definition_regex", "$^"))
    mdir = outdir / "markers"
    mdir.mkdir(exist_ok=True)
    handles = {mid: open(mdir / f"{mid}.faa", "w") for mid in matcher.ids}
    unmatched: Counter = Counter()
    n_rec = n_keep = n_noseq = 0
    rows = []
    try:
        for rec in parse_records(gbffs):
            n_rec += 1
            if n_rec % 2000 == 0:
                log.info("  parsed %d records (%d kept)", n_rec, n_keep)
            definition = rec.description
            length = len(rec.seq)
            try:
                gc = gc_fraction(rec.seq)
            except Exception:  # CON records without sequence content
                n_noseq += 1
                continue
            if filt.get("require_complete", True) and "complete genome" not in definition.lower():
                continue
            if excl.search(definition) or length < filt.get("min_length", 0) or length > filt.get("max_length", 1e12):
                continue
            organism = rec.annotations.get("organism", "")
            lineage = rec.annotations.get("taxonomy", [])
            grp = mapper.group_of(lineage, organism)
            found: dict[str, str] = {}
            n_cds = n_trna = n_rrna = 0
            tables: Counter = Counter()
            for f in rec.features:
                if f.type == "tRNA":
                    n_trna += 1
                elif f.type == "rRNA":
                    n_rrna += 1
                elif f.type == "CDS":
                    n_cds += 1
                    if "transl_table" in f.qualifiers:
                        tables[f.qualifiers["transl_table"][0]] += 1
                    gene = f.qualifiers.get("gene", [None])[0]
                    product = f.qualifiers.get("product", [None])[0]
                    mid = matcher.match(gene, product)
                    if mid is None:
                        if not matcher.is_ignored(product):
                            unmatched[(gene or "", product or "")] += 1
                        continue
                    aa = _protein(f, rec)
                    if aa and (mid not in found or len(aa) > len(found[mid])):
                        found[mid] = aa
            acc = rec.id
            for mid, aa in found.items():
                handles[mid].write(f">{acc}\n{aa}\n")
            rows.append({
                "accession": acc, "organism": organism, "taxid": _taxid(rec), "group": grp["id"],
                "supergroup": grp.get("supergroup", ""), "lineage": "; ".join(lineage), "definition": definition,
                "length": length, "gc": f"{gc:.4f}", "n_cds": n_cds, "n_trna": n_trna, "n_rrna": n_rrna,
                "transl_table": tables.most_common(1)[0][0] if tables else ("1" if n_cds else ""), "n_markers": len(found),
                "markers": ",".join(sorted(found)), "topology": rec.annotations.get("topology", ""),
                "date": rec.annotations.get("date", ""),
            })
            n_keep += 1
    finally:
        for h in handles.values():
            h.close()
    with open(genomes_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=GENOME_COLUMNS, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    with open(outdir / "unmatched_products.tsv", "w") as fh:
        fh.write("count\tgene\tproduct\n")
        for (g, p), c in unmatched.most_common():
            fh.write(f"{c}\t{g}\t{p}\n")
    marker_counts = {mid: sum(1 for r in rows if mid in r["markers"].split(",")) for mid in matcher.ids}
    write_manifest(outdir, ih, n_records=n_rec, n_kept=n_keep, n_no_sequence=n_noseq, marker_counts=marker_counts, gbff=[p.name for p in gbffs])
    log.info("extract-mito: %d records parsed, %d kept; markers: %s", n_rec, n_keep, marker_counts)
    return genomes_tsv


def load_genomes(cfg: Config) -> list[dict]:
    p = cfg.workdir / "mito" / "genomes.tsv"
    if not p.exists():
        raise SystemExit("run `mitophy extract-mito` first")
    with open(p) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["length"] = int(r["length"]); r["gc"] = float(r["gc"]); r["n_markers"] = int(r["n_markers"])
        r["n_cds"] = int(r["n_cds"]); r["n_trna"] = int(r["n_trna"]); r["n_rrna"] = int(r["n_rrna"])
        r["marker_set"] = set(r["markers"].split(",")) if r["markers"] else set()
    return rows

def read_marker_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    if not path.exists():
        return seqs
    with open(path) as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                seqs[name] = ""
            elif name:
                seqs[name] += line.strip()
    return seqs
