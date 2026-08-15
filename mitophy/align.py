"""Stage 6: per-marker alignment (MAFFT), trimming (trimAl), occupancy filtering and concatenation.

Inputs per analysis:
  origin          work/origin/markers/<m>.faa   (mito + bacterial orthologs, from `orthologs`)
  diversification work/mito/markers/<m>.faa filtered to work/diversification/selected.tsv
Outputs: work/<analysis>/aln/<m>.{faa,mafft.faa,trim.faa}, concat.faa, partitions.nex, occupancy.tsv, taxa.tsv
"""
from __future__ import annotations
import csv, json, logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from .config import Config
from .cache import hash_inputs, up_to_date, write_manifest
from .extract_mito import read_marker_fasta
from .tools import run as sh, which_first

log = logging.getLogger("mitophy")

def write_fasta(seqs: dict[str, str], path: Path, width: int = 80) -> None:
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), width):
                fh.write(v[i:i + width] + "\n")

def _marker_inputs(cfg: Config, analysis: str) -> tuple[dict[str, dict[str, str]], list[dict]]:
    """Return {marker: {taxon_label: seq}} and taxa metadata rows."""
    if analysis == "diversification":
        from .sampling import load_selected
        sel = load_selected(cfg)
        keep = {r["accession"] for r in sel}
        markers = {}
        for m in cfg.markers_for("diversification"):
            seqs = read_marker_fasta(cfg.workdir / "mito" / "markers" / f"{m['id']}.faa")
            markers[m["id"]] = {k: v for k, v in seqs.items() if k in keep}
        taxa = [{"label": r["accession"], "name": r["organism"], "group": r["group"], "supergroup": r["supergroup"],
                 "taxid": r["taxid"], "accession": r["accession"], "lineage": r["lineage"], "kind": "mito"} for r in sel]
    else:
        tdir = cfg.workdir / "origin"
        taxa_tsv = tdir / "taxa.tsv"
        if not taxa_tsv.exists():
            raise SystemExit("run `mitophy orthologs` first")
        with open(taxa_tsv) as fh:
            taxa = list(csv.DictReader(fh, delimiter="\t"))
        markers = {m["id"]: read_marker_fasta(tdir / "markers" / f"{m['id']}.faa") for m in cfg.markers_for("origin")}
    return markers, taxa

def _align_one(cfg: Config, analysis: str, mid: str, seqs: dict[str, str], adir: Path) -> Path | None:
    inp = adir / f"{mid}.faa"; aln = adir / f"{mid}.mafft.faa"; trim = adir / f"{mid}.trim.faa"
    if len(seqs) < 4:
        log.warning("%s/%s: only %d sequences, skipping", analysis, mid, len(seqs))
        return None
    write_fasta(seqs, inp)
    margs = cfg.tools["mafft"][analysis]
    sh(["mafft", *margs, "--thread", "1", inp], stdout_to=aln, log_file=adir / f"{mid}.mafft.log")
    targs = cfg.tools["trimal"][analysis]
    sh(["trimal", "-in", aln, "-out", trim, *targs], log_file=adir / f"{mid}.trimal.log")
    return trim

def run(cfg: Config, analysis: str) -> Path:
    outdir = cfg.sub(analysis)
    adir = outdir / "aln"; adir.mkdir(exist_ok=True)
    markers, taxa = _marker_inputs(cfg, analysis)
    occ = cfg.pipeline["occupancy"][analysis]
    ih = hash_inputs(objs=[{m: sorted(s.items()) for m, s in markers.items()}, cfg.tools["mafft"][analysis], cfg.tools["trimal"][analysis], occ])
    concat = outdir / "concat.faa"
    if up_to_date(outdir / "align", ih, [concat, outdir / "partitions.nex"], cfg.force):
        log.info("align/%s up to date", analysis); return concat
    threads = cfg.params.get("threads", 4)
    n_taxa_total = len({t for s in markers.values() for t in s})
    todo = {m: s for m, s in markers.items() if len(s) >= max(4, occ["min_taxa_per_marker_frac"] * n_taxa_total)}
    dropped_markers = sorted(set(markers) - set(todo))
    if dropped_markers:
        log.info("%s: dropping low-occupancy markers: %s", analysis, dropped_markers)
    with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
        futs = {m: ex.submit(_align_one, cfg, analysis, m, s, adir) for m, s in todo.items()}
        trimmed = {m: f.result() for m, f in futs.items()}
    trimmed = {m: p for m, p in trimmed.items() if p}
    alns = {m: read_marker_fasta(p) for m, p in trimmed.items()}
    # drop taxa with too few markers
    all_taxa = sorted({t for s in alns.values() for t in s})
    counts = {t: sum(1 for s in alns.values() if t in s and set(s[t]) - {"-", "X", "?"}) for t in all_taxa}
    kept_taxa = [t for t in all_taxa if counts[t] >= occ["min_markers_per_taxon"]]
    log.info("%s: %d markers, %d/%d taxa kept (min %d markers)", analysis, len(alns), len(kept_taxa), len(all_taxa), occ["min_markers_per_taxon"])
    parts, concat_seqs, pos = [], {t: [] for t in kept_taxa}, 1
    for m in sorted(alns):
        s = alns[m]
        L = len(next(iter(s.values())))
        for t in kept_taxa:
            concat_seqs[t].append(s.get(t, "-" * L))
        parts.append((m, pos, pos + L - 1)); pos += L
    write_fasta({t: "".join(v) for t, v in concat_seqs.items()}, concat)
    with open(outdir / "partitions.nex", "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for m, a, b in parts:
            fh.write(f"  charset {m} = {a}-{b};\n")
        fh.write("end;\n")
    with open(outdir / "occupancy.tsv", "w") as fh:
        fh.write("taxon\t" + "\t".join(m for m, _, _ in parts) + "\tn_markers\n")
        for t in kept_taxa:
            row = [("1" if (t in alns[m] and set(alns[m][t]) - {"-", "X", "?"}) else "0") for m, _, _ in parts]
            fh.write(f"{t}\t" + "\t".join(row) + f"\t{row.count('1')}\n")
    tmeta = [t for t in taxa if t["label"] in set(kept_taxa)]
    with open(outdir / "taxa.json", "w") as fh:
        json.dump(tmeta, fh, indent=1)
    write_manifest(outdir / "align", ih, n_taxa=len(kept_taxa), n_sites=pos - 1, markers=[m for m, _, _ in parts],
                   dropped_markers=dropped_markers, dropped_taxa=[t for t in all_taxa if t not in kept_taxa])
    return concat
