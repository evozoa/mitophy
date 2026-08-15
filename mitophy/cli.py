"""Command-line entry point: `mitophy <stage> [--profile P] [--force] [--workdir DIR]`."""
from __future__ import annotations
import argparse, logging, sys, time
from .config import load_config

STAGES = ["fetch-mito", "extract-mito", "fetch-bact", "orthologs", "sample", "align", "trees", "stats", "provenance", "site"]

def _run_stage(name: str, cfg):
    from . import fetch_mito, extract_mito, fetch_bact, orthologs, sampling, align, trees, stats, provenance
    from .site import build as site_build
    t0 = time.time()
    logging.getLogger("mitophy").info("=== stage %s (profile=%s) ===", name, cfg.profile)
    if name == "fetch-mito":
        fetch_mito.run(cfg)
    elif name == "extract-mito":
        extract_mito.run(cfg)
    elif name == "fetch-bact":
        fetch_bact.run(cfg)
    elif name == "orthologs":
        orthologs.run(cfg)
    elif name == "sample":
        sampling.run(cfg)
    elif name == "align":
        align.run(cfg, "origin"); align.run(cfg, "diversification")
    elif name == "trees":
        trees.run(cfg, "origin"); trees.run(cfg, "diversification")
    elif name == "stats":
        stats.run(cfg)
    elif name == "provenance":
        provenance.run(cfg)
    elif name == "site":
        site_build.run(cfg)
    else:
        raise SystemExit(f"unknown stage {name}")
    logging.getLogger("mitophy").info("=== stage %s done in %.0fs ===", name, time.time() - t0)

def main(argv=None):
    ap = argparse.ArgumentParser(prog="mitophy", description="Evolution of mitochondria: pipeline + site")
    ap.add_argument("stage", choices=STAGES + ["all", "status"])
    ap.add_argument("--profile", default="default")
    ap.add_argument("--workdir")
    ap.add_argument("--resultsdir")
    ap.add_argument("--force", action="store_true", help="ignore cached manifests and rerun")
    ap.add_argument("--threads", type=int)
    ap.add_argument("--analysis", choices=["origin", "diversification"], help="restrict align/trees to one analysis")
    ap.add_argument("--from-stage", dest="from_stage", choices=STAGES, help="with `all`: start from this stage")
    ap.add_argument("--to-stage", dest="to_stage", choices=STAGES, help="with `all`: stop after this stage")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S", stream=sys.stderr)
    cfg = load_config(a.profile, workdir=a.workdir, resultsdir=a.resultsdir, force=a.force, threads=a.threads)
    if a.stage == "status":
        from .cache import read_manifest
        for s in ["mito", "bact", "origin", "diversification"]:
            m = read_manifest(cfg.workdir / s)
            print(f"{s:16s} {m.get('finished', '-'):22s} {m.get('input_hash', '')[:12]}")
        return
    if a.stage == "all":
        stages = STAGES[:-1]  # everything except site (site is built separately in CI after committing results)
        if a.from_stage:
            stages = stages[stages.index(a.from_stage):]
        if a.to_stage:
            stages = stages[: stages.index(a.to_stage) + 1]
    else:
        stages = [a.stage]
    for s in stages:
        if s in ("align", "trees") and a.analysis:
            from . import align, trees
            (align if s == "align" else trees).run(cfg, a.analysis)
        else:
            _run_stage(s, cfg)

if __name__ == "__main__":
    main()
