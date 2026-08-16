"""Stage 4: assemble per-marker ortholog sets for the origin analysis (mito proteins + bacterial homologs).

1. mito: proteins for the taxa listed in config/taxa_origin.yml (name-based extraction, from `extract-mito`)
2. bacterial seeds: proteins whose RefSeq product name matches the marker's `bact.product` regexes
3. per-marker profile HMM built (pyhmmer) from an alignment of mito + seed sequences
4. hmmsearch against every proteome; best hit per proteome (bitscore/E-value thresholds), 2nd-best margin -> ambiguity flag;
   a protein claimed by two markers is kept for the marker where it scores higher (paralog resolution)
5. mito sequences are validated against the same HMM (drops mislabelled genes)
Outputs: work/origin/markers/<m>.faa, work/origin/hits.tsv, work/origin/taxa.tsv, work/origin/hmm/<m>.hmm
"""
from __future__ import annotations
import csv, json, logging, re
from collections import defaultdict
from pathlib import Path
import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock, MSAFile
from .config import Config
from .cache import hash_inputs, up_to_date, write_manifest
from .extract_mito import load_genomes, read_marker_fasta
from .fetch_bact import load_assemblies
from .align import write_fasta
from .tools import run as sh

log = logging.getLogger("mitophy")
ALPHA = Alphabet.amino()

def _s(x) -> str:
    return x.decode() if isinstance(x, bytes) else str(x)

def origin_mito_taxa(cfg: Config, genomes: list[dict]) -> list[dict]:
    by_name = defaultdict(list)
    for g in genomes:
        by_name[g["organism"]].append(g)
    out = []
    for e in cfg.taxa_origin["mito"]:
        cands = by_name.get(e["name"], [])
        if not cands and "accession" in e:
            cands = [g for g in genomes if g["accession"].split(".")[0] == e["accession"].split(".")[0]]
        if not cands:  # genus + epithet prefix match
            cands = [g for g in genomes if g["organism"].startswith(e["name"])]
        if not cands:
            log.warning("origin mito taxon not found in RefSeq: %s", e["name"]); continue
        best = max(cands, key=lambda g: (g["n_markers"], g["length"]))
        out.append(best)
    maxm = cfg.params["origin"].get("max_mito")
    if maxm:
        out = sorted(out, key=lambda g: -g["n_markers"])[:maxm]
    return out

def _read_proteome(path: Path, gcf: str) -> tuple[list[TextSequence], dict[str, str]]:
    seqs, products = [], {}
    name = None; buf = []
    def flush():
        if name is not None:
            seqs.append(TextSequence(name=f"{gcf}|{name}".encode(), sequence="".join(buf)))
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(); buf = []
                hdr = line[1:].strip()
                name = hdr.split()[0]
                desc = hdr[len(name):].strip()
                products[name] = re.sub(r"\s*\[[^\]]*\]$", "", desc)
            else:
                buf.append(line.strip())
        flush()
    return seqs, products

def run(cfg: Config) -> Path:
    outdir = cfg.sub("origin")
    genomes = load_genomes(cfg)
    mito = origin_mito_taxa(cfg, genomes)
    assemblies = load_assemblies(cfg)
    markers = cfg.markers_for("origin")
    ocfg = cfg.pipeline["orthologs"]
    pdir = cfg.workdir / "bact" / "proteomes"
    ih = hash_inputs(files=[pdir / f"{a['accession']}.faa" for a in assemblies] + [cfg.root / "config" / "markers.yml"],
                     objs=[[g["accession"] for g in mito], ocfg, [m["id"] for m in markers]])
    taxa_tsv = outdir / "taxa.tsv"
    if up_to_date(outdir / "orthologs", ih, [taxa_tsv], cfg.force):
        log.info("orthologs up to date"); return taxa_tsv
    mdir = outdir / "markers"; mdir.mkdir(exist_ok=True)
    hdir = outdir / "hmm"; hdir.mkdir(exist_ok=True)
    tmp = outdir / "tmp"; tmp.mkdir(exist_ok=True)
    threads = cfg.params.get("threads", 4)

    # ---- load proteomes -------------------------------------------------------------------------
    all_text: list[TextSequence] = []
    products: dict[str, dict[str, str]] = {}
    for a in assemblies:
        p = pdir / f"{a['accession']}.faa"
        if not p.exists():
            log.warning("missing proteome %s", p); continue
        seqs, prods = _read_proteome(p, a["accession"])
        all_text.extend(seqs); products[a["accession"]] = prods
    seq_by_name = {_s(s.name): s for s in all_text}
    block = DigitalSequenceBlock(ALPHA, [s.digitize(ALPHA) for s in all_text])
    log.info("orthologs: %d proteins from %d proteomes; %d mito genomes; %d markers", len(all_text), len(products), len(mito), len(markers))
    mito_acc = [g["accession"] for g in mito]

    # ---- build HMMs (seed alignments run with --thread 1: multithreaded MAFFT is not deterministic) ----------------
    builder = pyhmmer.plan7.Builder(ALPHA)
    background = pyhmmer.plan7.Background(ALPHA)
    hmms, mito_seqs = [], {}
    jobs = []
    for m in markers:
        mid = m["id"]
        ms = read_marker_fasta(cfg.workdir / "mito" / "markers" / f"{mid}.faa")
        ms = {k: v for k, v in ms.items() if k in mito_acc}
        mito_seqs[mid] = ms
        pats = [re.compile(p, re.I) for p in m.get("bact", {}).get("product", [])]
        med = sorted(len(v) for v in ms.values())[len(ms) // 2] if ms else 0
        seeds = {}
        for gcf, prods in products.items():
            for wp, prod in prods.items():
                if any(p.search(prod) for p in pats):
                    s = seq_by_name[f"{gcf}|{wp}"]
                    if len(s.sequence) >= ocfg["seed_min_len_frac"] * med:
                        seeds[f"{gcf}|{wp}"] = s.sequence
                        break  # one seed per proteome
        train = {**{f"m|{k}": v for k, v in ms.items()}, **{f"b|{k}": v for k, v in seeds.items()}}
        if len(train) < 3:
            log.warning("%s: too few training sequences (%d), skipping marker", mid, len(train)); continue
        fa = tmp / f"{mid}.seed.faa"; write_fasta(train, fa)
        jobs.append((mid, fa, tmp / f"{mid}.seed.aln", len(ms), len(seeds)))
    from concurrent.futures import ThreadPoolExecutor
    def _aln(job):
        mid, fa, aln, *_ = job
        sh(["mafft", "--auto", "--anysymbol", "--quiet", "--thread", "1", fa], stdout_to=aln, log_file=tmp / f"{mid}.mafft.log")
        return job
    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        done = list(ex.map(_aln, jobs))
    for mid, fa, aln, n_m, n_b in done:
        with MSAFile(str(aln), digital=True, alphabet=ALPHA) as f:
            msa = f.read()
        msa.name = mid.encode()
        hmm, _, _ = builder.build_msa(msa, background)
        with open(hdir / f"{mid}.hmm", "wb") as fh:
            hmm.write(fh)
        hmms.append(hmm)
        log.info("  %s: %d mito + %d bacterial seeds -> HMM (M=%d)", mid, n_m, n_b, hmm.M)

    # ---- search proteomes -----------------------------------------------------------------------
    hits_rows = []
    best: dict[str, dict[str, list]] = defaultdict(dict)  # marker -> gcf -> [(score, name, evalue), ...] sorted desc
    for top in pyhmmer.hmmsearch(hmms, block, cpus=threads, E=float(ocfg["evalue"])):
        mid = _s(top.query.name)
        per_gcf: dict[str, list] = defaultdict(list)
        for h in top:
            if h.score < ocfg["min_bitscore"]:
                continue
            gcf = _s(h.name).split("|")[0]
            per_gcf[gcf].append((h.score, _s(h.name), h.evalue))
        for gcf, lst in per_gcf.items():
            best[mid][gcf] = sorted(lst, reverse=True)
    # paralog resolution across markers: each protein assigned to the marker where it scores highest
    claims: dict[str, list[tuple[float, str, str]]] = defaultdict(list)  # protein -> [(score, marker, gcf)]
    chosen: dict[str, dict[str, tuple]] = defaultdict(dict)
    for mid, d in best.items():
        for gcf, lst in d.items():
            chosen[mid][gcf] = lst[0]
    changed = True
    while changed:
        changed = False
        claims.clear()
        for mid, d in chosen.items():
            for gcf, (score, name, ev) in d.items():
                claims[name].append((score, mid, gcf))
        for name, lst in claims.items():
            if len(lst) > 1:
                lst.sort(reverse=True)
                for score, mid, gcf in lst[1:]:
                    # give this marker its next-best hit in that proteome (if any)
                    cands = [c for c in best[mid][gcf] if c[1] != name and not any(c[1] == n2 for n2 in [])]
                    cands = [c for c in cands if all(chosen[m2].get(gcf, (0, None))[1] != c[1] for m2 in chosen)]
                    if cands:
                        chosen[mid][gcf] = cands[0]
                    else:
                        del chosen[mid][gcf]
                    changed = True
    # ---- validate mito sequences against HMMs & write marker FASTA -------------------------------
    n_out = {}
    for hmm in hmms:
        mid = _s(hmm.name)
        ms = mito_seqs[mid]
        mblock = DigitalSequenceBlock(ALPHA, [TextSequence(name=k.encode(), sequence=v).digitize(ALPHA) for k, v in ms.items()])
        ok = {_s(h.name) for h in next(iter(pyhmmer.hmmsearch([hmm], mblock, cpus=threads, E=1e-3)))} if ms else set()
        dropped = sorted(set(ms) - ok)
        if dropped:
            log.info("  %s: %d mito sequences failed HMM validation: %s", mid, len(dropped), ", ".join(dropped[:8]))
        out = {k: v for k, v in ms.items() if k in ok}
        for gcf, (score, name, ev) in chosen.get(mid, {}).items():
            out[gcf] = seq_by_name[name].sequence
            lst = best[mid][gcf]
            second = next((c for c in lst if c[1] != name), None)
            amb = bool(second and second[0] >= score * (1 - ocfg["ambiguity_margin"]))
            hits_rows.append({"marker": mid, "assembly": gcf, "protein": name.split("|", 1)[1], "bitscore": f"{score:.1f}",
                              "evalue": f"{ev:.1e}", "product": products[gcf].get(name.split("|", 1)[1], ""),
                              "second_best": f"{second[0]:.1f}" if second else "", "ambiguous": int(amb)})
        write_fasta(out, mdir / f"{mid}.faa")
        n_out[mid] = {"mito": len(ok), "bact": len(chosen.get(mid, {}))}
    with open(outdir / "hits.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["marker", "assembly", "protein", "bitscore", "evalue", "product", "second_best", "ambiguous"], delimiter="\t")
        w.writeheader(); w.writerows(sorted(hits_rows, key=lambda r: (r["marker"], r["assembly"])))
    # ---- taxa table ------------------------------------------------------------------------------
    with open(taxa_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["label", "name", "group", "supergroup", "kind", "accession", "taxid", "lineage"])
        for g in mito:
            w.writerow([g["accession"], g["organism"], g["group"], g["supergroup"], "mito", g["accession"], g["taxid"], g["lineage"]])
        for a in assemblies:
            w.writerow([a["accession"], a["name"], a["group"], a["supergroup"], "bacteria", a["accession"], "", a["config_group"]])
    write_manifest(outdir / "orthologs", ih, markers=n_out, n_mito=len(mito), n_bact=len(assemblies),
                   n_ambiguous=sum(r["ambiguous"] for r in hits_rows))
    log.info("orthologs: %s", n_out)
    return taxa_tsv
