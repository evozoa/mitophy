"""Stage 9: record run provenance (results/runs/<date>.json, results/latest.json) and prepend a CHANGELOG entry."""
from __future__ import annotations
import csv, json, logging, time
from pathlib import Path
from .config import Config
from .cache import read_manifest
from .tools import tool_versions
from . import __version__

log = logging.getLogger("mitophy")

def _sampled(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with open(path) as fh:
        return {r["accession"] for r in csv.DictReader(fh, delimiter="\t")}

def collect(cfg: Config) -> dict:
    wd, rd = cfg.workdir, cfg.resultsdir
    src = json.loads((wd / "mito" / "source.json").read_text()) if (wd / "mito" / "source.json").exists() else {}
    mito_m = read_manifest(wd / "mito")
    div_tree = json.loads((rd / "diversification" / "tree.json").read_text()) if (rd / "diversification" / "tree.json").exists() else {}
    ori_tree = json.loads((rd / "origin" / "tree.json").read_text()) if (rd / "origin" / "tree.json").exists() else {}
    summary = json.loads((rd / "diversification" / "summary.json").read_text()) if (rd / "diversification" / "summary.json").exists() else {}
    def treeinfo(t):
        return {k: t.get(k) for k in ("tool", "model", "lnL", "n_taxa", "n_sites", "markers", "rooted_on", "date", "runtime_s", "profile", "versions")}
    return {
        "date": time.strftime("%Y-%m-%d"), "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "profile": cfg.profile, "mitophy_version": __version__,
        "refseq_mito": {"files": {k: v.get("last_modified") for k, v in src.get("files", {}).items()},
                        "n_records": mito_m.get("n_records"), "n_complete_genomes": mito_m.get("n_kept")},
        "tools": tool_versions(cfg.tools),
        "origin": treeinfo(ori_tree) | {"orthologs": read_manifest(wd / "origin" / "orthologs").get("markers")},
        "diversification": treeinfo(div_tree) | {"n_sampled": summary.get("n_sampled")},
        "summary": summary,
    }

def run(cfg: Config) -> Path:
    rd = cfg.resultsdir; runs = rd / "runs"; runs.mkdir(parents=True, exist_ok=True)
    rec = collect(cfg)
    prev_files = sorted(runs.glob("*.json"))
    prev = json.loads(prev_files[-1].read_text()) if prev_files else None
    # sampled accession diff vs previous run
    cur = _sampled(rd / "diversification" / "genomes.tsv")
    prev_set = set(prev.get("sampled_accessions", [])) if prev else set()
    rec["sampled_accessions"] = sorted(cur)
    rec["diff"] = {"added": sorted(cur - prev_set), "removed": sorted(prev_set - cur)} if prev else {"added": sorted(cur), "removed": []}
    out = runs / f"{rec['date']}.json"
    out.write_text(json.dumps(rec, indent=1, default=str))
    (rd / "latest.json").write_text(json.dumps({k: v for k, v in rec.items() if k != "sampled_accessions"}, indent=1, default=str))
    # changelog
    cl = cfg.resultsdir.parent / "CHANGELOG.md"
    lines = [f"## {rec['date']} — pipeline run (profile `{cfg.profile}`)", ""]
    rs = rec["refseq_mito"]
    lines.append(f"- RefSeq mitochondrion release files: {', '.join(f'{k} ({v})' for k, v in rs['files'].items()) or 'n/a'}; "
                 f"{rs.get('n_complete_genomes')} complete mitogenomes parsed")
    o, d = rec["origin"], rec["diversification"]
    if o.get("tool"):
        lines.append(f"- Origin tree: {o['tool']} {o.get('model')}, {o.get('n_taxa')} taxa × {o.get('n_sites')} sites, lnL {o.get('lnL')}, rooted on {o.get('rooted_on')}")
    if d.get("tool"):
        lines.append(f"- Diversification tree: {d['tool']} {d.get('model')}, {d.get('n_taxa')} taxa × {d.get('n_sites')} sites, lnL {d.get('lnL')}")
    lines.append(f"- Sampled genomes: {len(cur)} (+{len(rec['diff']['added'])} / −{len(rec['diff']['removed'])} vs previous run)")
    lines.append(f"- Tools: " + ", ".join(f"{k} {v}" for k, v in rec["tools"].items()))
    entry = "\n".join(lines) + "\n\n"
    intro = "# Changelog\n\nAutomatically generated entries for pipeline runs; hand-written entries for content changes.\n\n"
    old = cl.read_text() if cl.exists() else intro
    i = old.find("\n## ")
    head, rest = (old[: i + 1], old[i + 1:]) if i >= 0 else (old if old.endswith("\n\n") else old + "\n\n", "")
    marker = f"## {rec['date']} — pipeline run"
    if marker in rest:
        import re
        rest = re.sub(rf"## {rec['date']} — pipeline run.*?(?=\n## |\Z)", entry.rstrip("\n") + "\n", rest, count=1, flags=re.S)
        cl.write_text(head + rest)
    else:
        cl.write_text(head + entry + rest)
    log.info("provenance written: %s", out)
    return out
