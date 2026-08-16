"""Small analyses computed at site-build time from the published trees/tables."""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import dendropy

def origin_headline(tree_json: dict) -> dict:
    """Where do the mitochondria branch? Returns monophyly, support and sister-group composition."""
    t = dendropy.Tree.get(data=tree_json["newick"], schema="newick", suppress_internal_node_taxa=True, preserve_underscores=True)
    t.is_rooted = True
    tips = tree_json["tips"]
    mito = [l for l, m in tips.items() if m.get("kind") == "mito"]
    if len(mito) < 2:
        return {"ok": False}
    mrca = t.mrca(taxon_labels=mito)
    clade = [n.taxon.label for n in mrca.leaf_iter()]
    mono = len(clade) == len(mito)
    def groups_of(labels):
        return Counter(tips.get(l, {}).get("group", "?") for l in labels)
    res = {"ok": True, "monophyletic": mono, "n_mito": len(mito), "clade_size": len(clade),
           "intruders": groups_of([l for l in clade if l not in set(mito)]) if not mono else {}}
    sup = tree_json.get("nodes", {}).get(mrca.label, {}) if mrca.label else {}
    res["support"] = sup
    if mrca.parent_node is not None:
        sisters = [c for c in mrca.parent_node.child_nodes() if c is not mrca]
        sis_labels = [n.taxon.label for s in sisters for n in s.leaf_iter()]
        res["sister_groups"] = groups_of(sis_labels).most_common()
        res["sister_size"] = len(sis_labels)
        res["sister_names"] = [tips[l]["name"] for l in sis_labels[:6]]
        psup = tree_json.get("nodes", {}).get(mrca.parent_node.label, {}) if mrca.parent_node.label else {}
        res["sister_support"] = psup
        # is the sister group exactly (all and only) Rickettsiales? or is mito sister to all alphaproteobacteria?
        sg = set(tips.get(l, {}).get("group") for l in sis_labels)
        alpha_groups = {g for l, m in tips.items() if m.get("supergroup") == "Alphaproteobacteria" for g in [m.get("group")]}
        n_alpha_total = sum(1 for m in tips.values() if m.get("supergroup") == "Alphaproteobacteria")
        n_alpha_in_sister = sum(1 for l in sis_labels if tips.get(l, {}).get("supergroup") == "Alphaproteobacteria")
        if sg == {"Rickettsiales"}:
            res["interpretation"] = "rickettsiales_sister"
        elif n_alpha_in_sister == n_alpha_total and n_alpha_total > 0 and all(tips.get(l, {}).get("supergroup") == "Alphaproteobacteria" for l in sis_labels):
            res["interpretation"] = "sister_to_alpha"
        elif "Rickettsiales" in sg:
            res["interpretation"] = "within_alpha_with_rickettsiales"
        else:
            res["interpretation"] = "other"
    return res

def occupancy_svg(occ_rows: list[dict], taxa: list[dict], out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "mitophy", "font.family": "sans-serif", "font.size": 8})
    tmeta = {t["label"]: t for t in taxa}
    markers = [k for k in occ_rows[0].keys() if k not in ("taxon", "n_markers")]
    rows = sorted(occ_rows, key=lambda r: (tmeta.get(r["taxon"], {}).get("kind", "") != "mito", tmeta.get(r["taxon"], {}).get("group", ""), r["taxon"]))
    import numpy as np
    mat = np.array([[int(r[m]) for m in markers] for r in rows])
    fig, ax = plt.subplots(figsize=(0.22 * len(markers) + 2.5, 0.16 * len(rows) + 1.2))
    ax.imshow(mat, cmap="Blues", vmin=0, vmax=1.4, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(markers))); ax.set_xticklabels(markers, rotation=90)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([(tmeta.get(r["taxon"], {}).get("name") or r["taxon"])[:32] for r in rows])
    for lab in ax.get_yticklabels():
        pass
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Marker occupancy (filled = present) — mitochondria first, then bacteria by group", loc="left", fontsize=9)
    fig.savefig(out, format="svg", bbox_inches="tight", transparent=True, metadata={"Date": None}); plt.close(fig)
