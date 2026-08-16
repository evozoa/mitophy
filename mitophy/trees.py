"""Stage 7: tree inference (IQ-TREE / FastTree), rooting, and tree.json export for the site."""
from __future__ import annotations
import json, logging, re, shutil, time
from pathlib import Path
import dendropy
from .config import Config
from .cache import hash_inputs, read_manifest, up_to_date, write_manifest
from .tools import run as sh, iqtree_bin, fasttree_bin, tool_versions
from .taxonomy import GroupMapper

log = logging.getLogger("mitophy")

OUTGROUP = {"origin": ["Betaproteobacteria", "Gammaproteobacteria", "Bacteria"], "diversification": ["Jakobida"]}

def _iqtree(cfg: Config, analysis: str, concat: Path, tdir: Path, p: dict) -> tuple[Path, dict]:
    ib = iqtree_bin(cfg.tools.get("iqtree", {}).get("bin", "auto"))
    threads = cfg.params.get("threads", 4)
    prefix = tdir / "iqtree"
    base = [ib, "-s", concat, "-T", str(threads), "--seed", str(cfg.pipeline.get("seed", 42)), "-redo", "--quiet"]
    model = p.get("model", "LG+F+G4")
    extra = list(p.get("iqtree_extra", []) or [])
    if p.get("pmsf_guide"):
        guide = tdir / "guide"
        sh([*base, "-m", p["pmsf_guide"], "--prefix", guide], log_file=tdir / "guide.stdout.log")
        cmd = [*base, "-m", model, "-ft", f"{guide}.treefile", "--prefix", prefix]
    else:
        cmd = [*base, "-m", model, "--prefix", prefix]
    if p.get("ufboot") and "--fast" in extra:
        log.warning("UFBoot is incompatible with --fast; disabling UFBoot for this run")
    elif p.get("ufboot"):
        cmd += ["-B", str(p["ufboot"])]
    if p.get("alrt"):
        cmd += ["-alrt", str(p["alrt"])]
    cmd += extra
    sh(cmd, log_file=tdir / "iqtree.stdout.log")
    info = {"tool": "IQ-TREE", "command": " ".join(map(str, cmd)), "model": model}
    rep = Path(f"{prefix}.iqtree")
    if rep.exists():
        txt = rep.read_text()
        m = re.search(r"Log-likelihood of the tree:\s*(-?[\d.]+)", txt); info["lnL"] = float(m.group(1)) if m else None
        m = re.search(r"Best-fit model according to \w+:\s*(\S+)", txt); info["model"] = m.group(1) if m else model
        m = re.search(r"Total wall-clock time used:\s*([\d.]+)", txt); info["wallclock_s"] = float(m.group(1)) if m else None
    return Path(f"{prefix}.treefile"), info

def _fasttree(cfg: Config, analysis: str, concat: Path, tdir: Path, p: dict) -> tuple[Path, dict]:
    fb = fasttree_bin(cfg.tools.get("fasttree", {}).get("bin", "auto"))
    out = tdir / "fasttree.nwk"
    args = p.get("fasttree_args", ["-lg", "-gamma"])
    import os, subprocess
    env = dict(os.environ, OMP_NUM_THREADS=str(cfg.params.get("threads", 4)))
    cmd = [fb, *args, str(concat)]
    log.info("$ %s", " ".join(cmd))
    with open(out, "w") as o, open(tdir / "fasttree.log", "w") as e:
        r = subprocess.run(cmd, stdout=o, stderr=e, env=env, text=True)
    if r.returncode != 0:
        raise RuntimeError("FastTree failed: " + (tdir / "fasttree.log").read_text()[-2000:])
    txt = (tdir / "fasttree.log").read_text()
    m = re.search(r"Gamma\(20\) LogLk\s*(-?[\d.]+)", txt) or re.search(r"LogLk\s*=\s*(-?[\d.]+)", txt)
    return out, {"tool": "FastTree", "command": " ".join(cmd), "model": " ".join(args).replace("-", "").upper() or "JTT",
                 "lnL": float(m.group(1)) if m else None}

def root_tree(newick_path: Path, taxa: list[dict], outgroup_groups: list[str]) -> tuple[dendropy.Tree, str]:
    tree = dendropy.Tree.get(path=str(newick_path), schema="newick", preserve_underscores=True, suppress_internal_node_taxa=True)
    tree.is_rooted = True
    og = [t["label"] for t in taxa if t.get("group") in outgroup_groups]
    labels = {n.taxon.label for n in tree.leaf_node_iter()}
    og = [l for l in og if l in labels]
    if not og:
        tree.reroot_at_midpoint(update_bipartitions=True, suppress_unifurcations=False)
        return tree, "midpoint (no outgroup present)"
    if len(og) == len(labels):
        tree.reroot_at_midpoint(update_bipartitions=True)
        return tree, "midpoint (all taxa are outgroup)"
    mrca = tree.mrca(taxon_labels=og)
    if mrca is tree.seed_node:  # outgroup paraphyletic w.r.t. current root: use complement
        ing = [l for l in labels if l not in set(og)]
        mrca = tree.mrca(taxon_labels=ing)
    if mrca is tree.seed_node or mrca.edge.length is None:
        tree.reroot_at_node(mrca, update_bipartitions=True) if mrca is not tree.seed_node else None
        how = "outgroup (paraphyletic; rooted at MRCA node)"
    else:
        tree.reroot_at_edge(mrca.edge, length1=mrca.edge.length / 2, length2=mrca.edge.length / 2, update_bipartitions=True)
        how = f"outgroup ({', '.join(outgroup_groups)})"
    return tree, how

def _support(label: str | None) -> dict:
    if not label:
        return {}
    parts = label.split("/")
    try:
        vals = [float(x) for x in parts]
    except ValueError:
        return {"label": label}
    if len(vals) == 2:
        return {"alrt": round(vals[0], 1), "ufboot": round(vals[1], 1)}
    v = vals[0]
    return {"support": round(v * 100 if v <= 1.0 else v, 1)}

def export_tree_json(tree: dendropy.Tree, taxa: list[dict], mapper: GroupMapper, info: dict, out: Path) -> dict:
    tmap = {t["label"]: t for t in taxa}
    # assign internal node ids and support
    nodes = {}
    for i, n in enumerate(tree.postorder_internal_node_iter()):
        nid = f"n{i}"
        sup = _support(n.label)
        n.label = nid
        if sup:
            nodes[nid] = sup
    tips = {}
    for n in tree.leaf_node_iter():
        lab = n.taxon.label
        t = tmap.get(lab, {"name": lab, "group": "Unassigned"})
        tips[lab] = {"name": t.get("name", lab), "group": t.get("group", ""), "supergroup": t.get("supergroup", ""),
                     "color": mapper.color(t.get("group", "")), "accession": t.get("accession", lab), "taxid": t.get("taxid", ""),
                     "kind": t.get("kind", "")}
    newick = tree.as_string(schema="newick", suppress_rooting=True, suppress_internal_node_labels=False,
                            unquoted_underscores=True, suppress_annotations=True).strip()
    groups = sorted({t["group"] for t in tips.values()})
    data = {"newick": newick, "tips": tips, "nodes": nodes, "legend": mapper.legend(groups), "n_taxa": len(tips), **info}
    out.write_text(json.dumps(data, indent=None, separators=(",", ":")))
    return data

def run(cfg: Config, analysis: str) -> Path:
    outdir = cfg.sub(analysis)
    concat = outdir / "concat.faa"
    if not concat.exists():
        raise SystemExit(f"run `mitophy align --analysis {analysis}` first")
    p = cfg.params[analysis]
    versions = tool_versions(cfg.tools)
    ih = hash_inputs(files=[concat], objs=[p, versions.get("iqtree"), versions.get("fasttree")])
    tdir = outdir / "tree"; tdir.mkdir(exist_ok=True)
    final_nwk, final_json = outdir / "tree.nwk", outdir / "tree.json"
    if up_to_date(tdir, ih, [final_nwk, final_json], cfg.force):
        log.info("trees/%s up to date", analysis); return final_nwk
    t0 = time.time()
    method = p.get("method", "iqtree")
    raw, info = (_iqtree if method == "iqtree" else _fasttree)(cfg, analysis, concat, tdir, p)
    taxa = json.loads((outdir / "taxa.json").read_text())
    tree, how = root_tree(raw, taxa, OUTGROUP[analysis])
    am = read_manifest(outdir / "align")
    info["support_type"] = ("SH-aLRT / UFBoot" if p.get("alrt") and p.get("ufboot") else "UFBoot" if p.get("ufboot") else "SH-aLRT" if p.get("alrt") else "none") if method == "iqtree" else "SH-like local support (FastTree)"
    info.update({"rooted_on": how, "n_sites": am.get("n_sites"), "markers": am.get("markers"), "analysis": analysis,
                 "date": time.strftime("%Y-%m-%d"), "runtime_s": round(time.time() - t0), "profile": cfg.profile,
                 "versions": {k: versions[k] for k in ("iqtree", "fasttree", "mafft", "trimal") if k in versions}})
    tree.write(path=str(final_nwk), schema="newick", suppress_rooting=True, unquoted_underscores=True)
    export_tree_json(tree, taxa, GroupMapper(cfg.supergroups), info, final_json)
    write_manifest(tdir, ih, **{k: v for k, v in info.items() if k != "markers"})
    # publish to results/
    rdir = cfg.resultsdir / (analysis + ("_deep" if cfg.profile == "deep" else "")); rdir.mkdir(parents=True, exist_ok=True)
    for f in [final_nwk, final_json, outdir / "occupancy.tsv", outdir / "partitions.nex", outdir / "taxa.json"]:
        shutil.copy(f, rdir / f.name)
    import gzip
    with open(concat, "rb") as src, gzip.GzipFile(rdir / "concat.faa.gz", "wb", mtime=0) as dst:   # mtime=0: reproducible bytes
        shutil.copyfileobj(src, dst)
    logf = tdir / ("iqtree.log" if method == "iqtree" else "fasttree.log")
    if logf.exists():
        with open(logf, "rb") as src, gzip.GzipFile(rdir / (logf.name + ".gz"), "wb", mtime=0) as dst:
            shutil.copyfileobj(src, dst)
    if (tdir / "iqtree.iqtree").exists():
        shutil.copy(tdir / "iqtree.iqtree", rdir / "iqtree.report.txt")
    (rdir / "manifest.json").write_text(json.dumps({**read_manifest(tdir), "align": am}, indent=2, default=str))
    log.info("trees/%s: %s (%s), lnL=%s, %s", analysis, info["tool"], info.get("model"), info.get("lnL"), how)
    return final_nwk
