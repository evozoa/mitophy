"""Stage 8: summary statistics of RefSeq mitochondrial genomes (all + sampled) and SVG charts for the site."""
from __future__ import annotations
import csv, json, logging, math
from collections import Counter, defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from .config import Config
from .cache import hash_inputs, up_to_date, write_manifest
from .extract_mito import load_genomes
from .sampling import load_selected

log = logging.getLogger("mitophy")
BLUE, INK, MUTED, GRID = "#2a78d6", "#0b0b0b", "#52514e", "#e5e4e0"
rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "mitophy", "font.family": "sans-serif", "font.size": 10, "axes.edgecolor": MUTED,
                 "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False,
                 "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                 "figure.facecolor": "none", "axes.facecolor": "none"})

def _save(fig, path: Path):
    fig.savefig(path, format="svg", bbox_inches="tight", transparent=True, dpi=200, metadata={"Date": None}); plt.close(fig)

def chart_length_by_group(rows: list[dict], out: Path, min_n: int = 5):
    by = defaultdict(list)
    for r in rows:
        by[r["group"]].append(r["length"])
    groups = [g for g in by if len(by[g]) >= min_n] or list(by)
    groups.sort(key=lambda g: (sorted(by[g])[len(by[g]) // 2], g))
    fig, ax = plt.subplots(figsize=(8, 0.32 * len(groups) + 1.2))
    import random
    rnd = random.Random(1)
    for i, g in enumerate(groups):
        ys = [i + rnd.uniform(-0.25, 0.25) for _ in by[g]]
        ax.scatter(by[g], ys, s=9, color=BLUE, alpha=0.45, linewidths=0, rasterized=True)
        med = sorted(by[g])[len(by[g]) // 2]
        ax.plot([med, med], [i - 0.35, i + 0.35], color=INK, lw=1.5)
    ax.set_xscale("log"); ax.set_yticks(range(len(groups))); ax.set_yticklabels([f"{g} (n={len(by[g])})" for g in groups])
    ax.set_xlabel("mitochondrial genome length (bp, log scale) — bar = median"); ax.grid(axis="y", visible=False)
    _save(fig, out)

def chart_genes_vs_length(rows: list[dict], out: Path, label_names: list[str]):
    fig, ax = plt.subplots(figsize=(7.5, 5))
    xs = [r["length"] for r in rows]; ys = [r["n_cds"] for r in rows]
    ax.scatter(xs, ys, s=12, color=BLUE, alpha=0.5, linewidths=0, rasterized=True)
    ax.set_xscale("log"); ax.set_xlabel("genome length (bp, log scale)"); ax.set_ylabel("annotated protein-coding genes")
    for r in rows:
        if r["organism"] in label_names:
            ax.annotate(r["organism"], (r["length"], r["n_cds"]), fontsize=8, color=INK, xytext=(4, 3), textcoords="offset points")
    _save(fig, out)

def chart_transl_tables(rows: list[dict], out: Path):
    c = Counter((r["supergroup"] or "Other", r["transl_table"] or "?") for r in rows)
    sgs = sorted({k[0] for k in c}, key=lambda s: (-sum(v for k, v in c.items() if k[0] == s), s))
    tables = sorted({k[1] for k in c}, key=lambda t: (t == "?", int(t) if t.isdigit() else 99))
    fig, ax = plt.subplots(figsize=(7.5, 0.5 * len(sgs) + 1.5))
    import numpy as np
    mat = np.array([[c.get((s, t), 0) for t in tables] for s in sgs], dtype=float)
    frac = mat / mat.sum(axis=1, keepdims=True)
    im = ax.imshow(frac, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(tables))); ax.set_xticklabels([f"table {t}" for t in tables], rotation=45, ha="right")
    ax.set_yticks(range(len(sgs))); ax.set_yticklabels([f"{s} (n={int(mat[i].sum())})" for i, s in enumerate(sgs)])
    for i in range(len(sgs)):
        for j in range(len(tables)):
            if mat[i, j] > 0:
                ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=7, color="#111111" if frac[i, j] < 0.6 else "white")
    ax.grid(False); ax.set_title("Genetic code (/transl_table) by supergroup — fraction of genomes", fontsize=10, loc="left")
    _save(fig, out)

def chart_gene_retention(rows: list[dict], markers: list[str], out: Path, min_n: int = 5):
    import numpy as np
    by = defaultdict(list)
    for r in rows:
        by[r["group"]].append(r["marker_set"])
    groups = [g for g in by if len(by[g]) >= min_n] or list(by)
    mat = np.array([[sum(1 for s in by[g] if m in s) / len(by[g]) for m in markers] for g in groups])
    order = sorted(range(len(groups)), key=lambda i: (-mat[i].sum(), groups[i]))
    groups = [groups[i] for i in order]; mat = mat[order]
    fig, ax = plt.subplots(figsize=(0.32 * len(markers) + 3, 0.32 * len(groups) + 1.5))
    ax.imshow(mat, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(markers))); ax.set_xticklabels(markers, rotation=90, fontsize=8)
    ax.set_yticks(range(len(groups))); ax.set_yticklabels([f"{g} (n={len(by[g])})" for g in groups], fontsize=8)
    ax.grid(False); ax.set_title("Fraction of genomes encoding each marker gene (RefSeq, all complete mitogenomes)", fontsize=10, loc="left")
    _save(fig, out)

def run(cfg: Config) -> Path:
    outdir = cfg.resultsdir / "diversification"; outdir.mkdir(parents=True, exist_ok=True)
    cdir = outdir / "charts"; cdir.mkdir(exist_ok=True)
    gen_tsv = cfg.workdir / "mito" / "genomes.tsv"; sel_tsv = cfg.workdir / "diversification" / "selected.tsv"
    ih = hash_inputs(files=[gen_tsv, sel_tsv, Path(__file__)])
    if up_to_date(cfg.sub("diversification") / "stats", ih, [outdir / "genomes.tsv", outdir / "stats_by_group.tsv"], cfg.force):
        log.info("stats up to date"); return outdir
    genomes = load_genomes(cfg); sel = load_selected(cfg)
    markers = [m["id"] for m in cfg.marker_list]
    # sampled genome table (published)
    cols = ["accession", "organism", "taxid", "group", "supergroup", "stratum", "length", "gc", "n_cds", "n_trna", "n_rrna", "transl_table", "n_markers", "markers"]
    with open(outdir / "genomes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(sel)
    # per-group aggregate over ALL RefSeq genomes
    by = defaultdict(list)
    for r in genomes:
        by[(r["supergroup"], r["group"])].append(r)
    with open(outdir / "stats_by_group.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["supergroup", "group", "n_genomes", "n_sampled", "median_length", "min_length", "max_length", "median_gc", "median_n_cds", "common_transl_table"] + [f"frac_{m}" for m in markers])
        nsel = Counter(r["group"] for r in sel)
        for (sg, g), rs in sorted(by.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            L = sorted(r["length"] for r in rs); G = sorted(r["gc"] for r in rs); C = sorted(r["n_cds"] for r in rs)
            tt = Counter(r["transl_table"] for r in rs).most_common(1)[0][0]
            w.writerow([sg, g, len(rs), nsel.get(g, 0), L[len(L) // 2], L[0], L[-1], f"{G[len(G) // 2]:.3f}", C[len(C) // 2], tt] +
                       [f"{sum(1 for r in rs if m in r['marker_set']) / len(rs):.3f}" for m in markers])
    summary = {"n_refseq_genomes": len(genomes), "n_sampled": len(sel), "n_groups": len(by),
               "n_by_supergroup": dict(Counter(r["supergroup"] for r in genomes).most_common()),
               "length_min": min(r["length"] for r in genomes), "length_max": max(r["length"] for r in genomes),
               "largest": max(genomes, key=lambda r: r["length"])["organism"], "smallest": min(genomes, key=lambda r: r["length"])["organism"],
               "most_genes": max(genomes, key=lambda r: r["n_cds"])["organism"], "most_genes_n": max(r["n_cds"] for r in genomes),
               "transl_tables": dict(Counter(r["transl_table"] for r in genomes).most_common())}
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    labels = ["Andalucia godoyi", "Reclinomonas americana", "Homo sapiens", "Marchantia polymorpha", "Plasmodium falciparum",
              "Saccharomyces cerevisiae", "Chondrus crispus", "Cucumis melo", "Silene conica", "Tetrahymena thermophila", "Chlamydomonas reinhardtii"]
    chart_length_by_group(genomes, cdir / "length_by_group.svg")
    chart_genes_vs_length(genomes, cdir / "genes_vs_length.svg", labels)
    chart_transl_tables(genomes, cdir / "transl_tables.svg")
    chart_gene_retention(genomes, [m for m in markers if m in {"cox1","cox2","cox3","cob","atp6","atp8","atp9","atp1","nad1","nad2","nad3","nad4","nad4L","nad5","nad6","nad7","nad8","nad9","nad11","sdh2","sdh3","sdh4","tufA","rpoB","secY","rps3","rps12","rpl2","rpl5","rpl16","rpl14","rps19","rps7","rps11","rpl6"}], cdir / "gene_retention.svg")
    write_manifest(cfg.sub("diversification") / "stats", ih, **summary)
    return outdir
