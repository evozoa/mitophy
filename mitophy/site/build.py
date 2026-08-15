"""Stage 10: render the static site into _site/ from results/, data/literature/ and config/."""
from __future__ import annotations
import csv, json, logging, shutil, time
from pathlib import Path
import yaml
from jinja2 import Environment, FileSystemLoader, select_autoescape
from ..config import Config, DATA_DIR
from .literature import load_literature, ReferenceFormatter
from .analysis import origin_headline, occupancy_svg

log = logging.getLogger("mitophy")
TEMPLATES = Path(__file__).parent / "templates"


def _read_json(p: Path, default=None):
    return json.loads(p.read_text()) if p.exists() else default


def _read_tsv(p: Path) -> list[dict]:
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_context(cfg: Config, out: Path) -> dict:
    rd = cfg.resultsdir
    refs = yaml.safe_load((DATA_DIR / "literature" / "references.yml").read_text())
    fmt = ReferenceFormatter(refs)
    lit = load_literature(DATA_DIR / "literature", out / "data" / "literature")
    hyps = yaml.safe_load((DATA_DIR / "literature" / "hypotheses.yml").read_text())
    latest = _read_json(rd / "latest.json", {})
    runs = sorted((_read_json(p) for p in (rd / "runs").glob("*.json")), key=lambda r: r["date"], reverse=True) if (rd / "runs").exists() else []
    origin = {"tree": _read_json(rd / "origin" / "tree.json"), "manifest": _read_json(rd / "origin" / "manifest.json", {}),
              "occupancy": _read_tsv(rd / "origin" / "occupancy.tsv"), "taxa": _read_json(rd / "origin" / "taxa.json", [])}
    origin_deep = {"tree": _read_json(rd / "origin_deep" / "tree.json"), "manifest": _read_json(rd / "origin_deep" / "manifest.json", {})}
    div = {"tree": _read_json(rd / "diversification" / "tree.json"), "manifest": _read_json(rd / "diversification" / "manifest.json", {}),
           "genomes": _read_tsv(rd / "diversification" / "genomes.tsv"), "by_group": _read_tsv(rd / "diversification" / "stats_by_group.tsv"),
           "summary": _read_json(rd / "diversification" / "summary.json", {}),
           "charts": sorted(p.name for p in (rd / "diversification" / "charts").glob("*.svg")) if (rd / "diversification" / "charts").exists() else []}
    # headline: where do mitochondria branch in the de novo origin tree?
    origin["headline"] = origin_headline(origin["tree"]) if origin["tree"] else None
    origin_deep["headline"] = origin_headline(origin_deep["tree"]) if origin_deep["tree"] else None
    if origin["occupancy"]:
        occupancy_svg(origin["occupancy"], origin["taxa"], out / "data" / "origin_occupancy.svg")
    # gene-content matrix for representative genomes (from the sampled genome table)
    reps = [n for n in cfg.sampling.get("must_include", [])]
    gen_by_name = {g["organism"]: g for g in div["genomes"]}
    matrix_markers = [m["id"] for m in cfg.marker_list]
    matrix_rows = []
    for n in reps:
        g = gen_by_name.get(n)
        if g:
            ms = set(g["markers"].split(",")) if g["markers"] else set()
            matrix_rows.append({"organism": n, "group": g["group"], "accession": g["accession"], "length": int(g["length"]),
                                "n_cds": g["n_cds"], "transl_table": g["transl_table"], "present": [m in ms for m in matrix_markers]})
    matrix_rows.sort(key=lambda r: -sum(r["present"]))
    # tables of taxa used in the origin analysis
    bact_rows = []
    for grp, entries in cfg.taxa_origin["bacteria"].items():
        for e in entries:
            bact_rows.append({**e, "config_group": grp})
    # inline SVG charts with theme-aware ink (matplotlib writes fixed hex colours; swap them for currentColor)
    def inline_svg(path: Path) -> str:
        if not path.exists():
            return ""
        svg = path.read_text()
        for hexc in ("#0b0b0b", "#52514e", "#222222", "#000000"):
            svg = svg.replace(hexc, "currentColor")
        svg = svg.replace("<svg ", '<svg style="max-width:100%;height:auto" ', 1)
        return svg
    charts_inline = {name: inline_svg(rd / "diversification" / "charts" / name) for name in div["charts"]}
    occupancy_inline = inline_svg(out / "data" / "origin_occupancy.svg") if origin["occupancy"] else ""
    ctx = {
        "charts_inline": charts_inline, "occupancy_inline": occupancy_inline,
        "site": cfg.site, "now": time.strftime("%Y-%m-%d"), "latest": latest, "runs": runs,
        "origin": origin, "origin_deep": origin_deep, "div": div, "lit": lit, "hyps": hyps, "refs": refs, "cite": fmt.cite, "ref_html": fmt.html,
        "markers": cfg.marker_list, "markers_origin": cfg.markers_for("origin"), "markers_div": cfg.markers_for("diversification"),
        "matrix_markers": matrix_markers, "matrix_rows": matrix_rows, "mito_taxa": cfg.taxa_origin["mito"], "bact_rows": bact_rows,
        "sampling": cfg.sampling, "pipeline": cfg.pipeline, "supergroups": cfg.supergroups["groups"],
        "transl_names": {"1": "standard", "2": "vertebrate mito", "3": "yeast mito", "4": "mold/protozoan/coelenterate mito (UGA=Trp)", "5": "invertebrate mito",
                         "9": "echinoderm/flatworm mito", "11": "bacterial/plastid", "13": "ascidian mito", "14": "alternative flatworm mito", "16": "chlorophycean mito",
                         "21": "trematode mito", "22": "Scenedesmus mito", "23": "Thraustochytrium mito", "24": "Rhabdopleuridae mito", "25": "SR1/Gracilibacteria", "33": "Cephalodiscidae mito"},
    }
    return ctx


def run(cfg: Config, out: Path | None = None) -> Path:
    out = Path(out) if out else cfg.root / "_site"
    if out.exists():
        shutil.rmtree(out)
    (out / "data").mkdir(parents=True)
    # static assets
    shutil.copytree(cfg.root / "static", out / "static")
    # results data for the viewer + downloads
    rd = cfg.resultsdir
    for sub in ("origin", "origin_deep", "diversification"):
        if (rd / sub).exists():
            shutil.copytree(rd / sub, out / "data" / sub)
    if (rd / "latest.json").exists():
        shutil.copy(rd / "latest.json", out / "data" / "latest.json")
    if (rd / "runs").exists():
        shutil.copytree(rd / "runs", out / "data" / "runs")
    ctx = build_context(cfg, out)
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=select_autoescape(["html"]), trim_blocks=True, lstrip_blocks=True)
    env.filters["fmtnum"] = lambda v: f"{int(v):,}" if v not in (None, "") and str(v).lstrip("-").isdigit() else v
    env.filters["fmtfloat"] = lambda v, d=1: (f"{float(v):,.{d}f}" if v not in (None, "") else "—")
    env.filters["tojson_compact"] = lambda v: json.dumps(v, separators=(",", ":"))
    def fmtsupport(v):
        if not v:
            return "n/a"
        parts = []
        if "alrt" in v: parts.append(f"SH-aLRT {v['alrt']:g}")
        if "ufboot" in v: parts.append(f"UFBoot {v['ufboot']:g}")
        if "support" in v: parts.append(f"{v['support']:g}")
        return " / ".join(parts) or json.dumps(v)
    env.filters["fmtsupport"] = fmtsupport
    for page in ("index", "origin", "diversification", "methods", "changelog"):
        tpl = env.get_template(f"{page}.html")
        (out / f"{page}.html").write_text(tpl.render(page=page, **ctx))
    (out / ".nojekyll").write_text("")
    log.info("site built in %s (%d files)", out, sum(1 for _ in out.rglob("*") if _.is_file()))
    return out
