"""Stage 5: deterministic, sticky, taxonomically stratified subsampling of mitochondrial genomes for the diversification tree."""
from __future__ import annotations
import csv, logging, re
from pathlib import Path
from .config import Config
from .cache import hash_inputs, up_to_date, write_manifest
from .extract_mito import load_genomes, GENOME_COLUMNS
from .taxonomy import lineage_list

log = logging.getLogger("mitophy")
_acc_num = re.compile(r"(\d+)")

def _acc_key(acc: str) -> tuple:
    m = _acc_num.search(acc)
    return (acc[:3], int(m.group(1)) if m else 0)

def _stratum(lineage: list[str], strata: list[dict]):
    for s in strata:
        if s["lineage_contains"] in lineage:
            return s
    return None

class _Node:
    __slots__ = ("children", "items")
    def __init__(self):
        self.children: dict[str, _Node] = {}
        self.items: list[dict] = []

def _build_trie(items: list[dict], root_tag: str) -> _Node:
    root = _Node()
    for g in items:
        lin = lineage_list(g["lineage"])
        i = lin.index(root_tag) if root_tag in lin else 0
        node = root
        for part in lin[i + 1:]:
            node = node.children.setdefault(part, _Node())
        node.items.append(g)
    return root

def _size(node: _Node) -> int:
    return len(node.items) + sum(_size(c) for c in node.children.values())

def balanced_pick(node: _Node, quota: int, score) -> list[dict]:
    """Pick up to `quota` items from the subtree, dividing the quota as evenly as possible among children
    (leftover from small children is redistributed to larger ones), preferring best `score` at the leaves."""
    if quota <= 0:
        return []
    avail = _size(node)
    if avail <= quota:
        out = list(node.items)
        for c in node.children.values():
            out.extend(balanced_pick(c, _size(c), score))
        return out
    # branches = direct items (each its own pseudo-branch) + children
    branches = [("item", g) for g in sorted(node.items, key=score)] + [("child", c) for c in node.children.values()]
    caps = [1 if k == "item" else _size(c) for k, c in branches]
    alloc = [0] * len(branches)
    remaining = quota
    active = [i for i in range(len(branches)) if caps[i] > 0]
    while remaining > 0 and active:
        share = max(1, remaining // len(active))
        progressed = False
        # small branches first so they are filled exactly; when only leftovers remain, prefer the largest branches.
        # Ties always keep the original branch order (items sorted by score first), so results are a fixed point
        # under sticky re-runs (a plain reverse sort would flip tie order and make the selection flip-flop).
        leftover = remaining < len(active) * share
        for i in sorted(active, key=lambda i: ((-(caps[i] - alloc[i]) if leftover else (caps[i] - alloc[i])), i)):
            give = min(share, caps[i] - alloc[i], remaining)
            if give > 0:
                alloc[i] += give; remaining -= give; progressed = True
            if remaining == 0:
                break
        active = [i for i in active if alloc[i] < caps[i]]
        if not progressed:
            break
    out = []
    for (kind, obj), n in zip(branches, alloc):
        if n <= 0:
            continue
        out.append(obj) if kind == "item" else out.extend(balanced_pick(obj, n, score))
    return out

def previous_selection(cfg: Config) -> set[str]:
    p = cfg.resultsdir / "diversification" / "genomes.tsv"
    if not p.exists():
        return set()
    with open(p) as fh:
        return {r["accession"] for r in csv.DictReader(fh, delimiter="\t")}

def select(genomes: list[dict], cfg: Config, previous: set[str] | None = None) -> list[dict]:
    scfg = cfg.sampling
    strata = scfg["strata"]
    min_markers = cfg.pipeline["occupancy"]["diversification"]["min_markers_per_taxon"]
    div_markers = {m["id"] for m in cfg.markers_for("diversification")}
    must = set(scfg.get("must_include", []))
    previous = previous if (previous is not None and scfg.get("sticky", True)) else set()
    max_taxa = cfg.params["diversification"].get("max_taxa")

    def score(g):  # lower is better
        return (0 if g["accession"] in previous else 1,
                -len(g["marker_set"] & div_markers),
                -g["length"], _acc_key(g["accession"]))

    # one genome per organism name (prefer best-scoring), to avoid many strains of one species
    by_org: dict[str, dict] = {}
    for g in genomes:
        if len(g["marker_set"] & div_markers) < min_markers:
            continue
        cur = by_org.get(g["organism"])
        if cur is None or score(g) < score(cur):
            by_org[g["organism"]] = g
    eligible = list(by_org.values())

    chosen: dict[str, dict] = {}
    for g in eligible:
        if g["organism"] in must or g["accession"] in must or g["accession"].split(".")[0] in must:
            g["stratum"] = "must_include"; g["group_key"] = g["organism"]
            chosen[g["accession"]] = g

    per_stratum: dict[str, list[dict]] = {}
    for g in eligible:
        if g["accession"] in chosen:
            continue
        st = _stratum(lineage_list(g["lineage"]), strata)
        if st is None:
            continue
        g["stratum"] = st["name"]
        per_stratum.setdefault(st["name"], []).append(g)

    picked_by_stratum: dict[str, list[dict]] = {}
    for s in strata:
        items = per_stratum.get(s["name"], [])
        trie = _build_trie(items, s["lineage_contains"])
        picked = balanced_pick(trie, s.get("total_cap", len(items)), score)
        # order picked so that trimming (pop from end) removes lowest-priority first: keep balanced order but sort by score within
        picked_by_stratum[s["name"]] = sorted(picked, key=score)
    for lst in picked_by_stratum.values():
        for g in lst:
            g["group_key"] = ";".join(lineage_list(g["lineage"])[-3:-1])
    total = len(chosen) + sum(len(v) for v in picked_by_stratum.values())
    if max_taxa and total > max_taxa:
        # trim round-robin from strata (removing lowest-priority = last picked) until within max_taxa
        order = sorted(picked_by_stratum, key=lambda k: -len(picked_by_stratum[k]))
        while total > max_taxa:
            for k in order:
                if picked_by_stratum[k] and total > max_taxa:
                    picked_by_stratum[k].pop(); total -= 1
            # (note: trimming breaks perfect balance slightly; prefer setting total_cap so that this rarely triggers)
    for lst in picked_by_stratum.values():
        for g in lst:
            chosen[g["accession"]] = g
    sel = sorted(chosen.values(), key=lambda g: (g["stratum"], g["group_key"], g["organism"]))
    log.info("sampling: %d eligible, %d selected (%s)", len(eligible), len(sel),
             ", ".join(f"{k}={len(v)}" for k, v in picked_by_stratum.items()) + f", must_include={sum(1 for g in sel if g['stratum']=='must_include')}")
    return sel

def run(cfg: Config) -> Path:
    outdir = cfg.sub("diversification")
    gen_tsv = cfg.workdir / "mito" / "genomes.tsv"
    ih = hash_inputs(files=[gen_tsv, cfg.root / "config" / "sampling_diversification.yml"],
                     objs=[cfg.params["diversification"], cfg.pipeline["occupancy"]["diversification"], sorted(previous_selection(cfg))])
    out = outdir / "selected.tsv"
    if up_to_date(outdir / "sample", ih, [out], cfg.force):
        log.info("sample up to date"); return out
    genomes = load_genomes(cfg)
    sel = select(genomes, cfg, previous_selection(cfg))
    cols = GENOME_COLUMNS + ["stratum", "group_key"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(sel)
    write_manifest(outdir / "sample", ih, n_selected=len(sel))
    return out

def load_selected(cfg: Config) -> list[dict]:
    p = cfg.workdir / "diversification" / "selected.tsv"
    if not p.exists():
        raise SystemExit("run `mitophy sample` first")
    with open(p) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["marker_set"] = set(r["markers"].split(",")) if r["markers"] else set()
    return rows
