import json
from pathlib import Path
import dendropy
from mitophy.config import load_config
from mitophy.extract_mito import MarkerMatcher, normalise
from mitophy.sampling import balanced_pick, _build_trie
from mitophy.trees import _support, root_tree
from mitophy.site.literature import load_literature

def test_config_profiles():
    c = load_config("quick")
    assert c.params["origin"]["max_mito"] == 12 and c.params["origin"]["ufboot"] == 0
    assert "--fast" in c.params["origin"]["iqtree_extra"]
    d = load_config()
    assert d.params["origin"]["ufboot"] == 1000

def test_marker_matching():
    m = MarkerMatcher(load_config().markers)
    assert m.match("COX1", None) == "cox1"
    assert m.match(None, "cytochrome c oxidase subunit I") == "cox1"
    assert m.match("ND4L", "NADH dehydrogenase subunit 4L") == "nad4L"
    assert m.match("MT-CYB", "cytochrome b") == "cob"
    assert m.match("nad9", None) == "nad9"
    assert m.match(None, "hypothetical protein") is None
    assert m.match(None, "ATP synthase F0 subunit 6") == "atp6"
    assert m.match("ymf16", "SecY-independent transporter protein") is None
    assert normalise("  Cytochrome-B ") == "cytochrome b"

def test_balanced_pick_is_balanced_and_deterministic():
    items = []
    for i in range(50):  # 50 vertebrates
        items.append({"accession": f"NC_{i:06d}", "lineage": "Eukaryota; Metazoa; Chordata; Vertebrata; Mammalia; Sp%d" % i, "n": i})
    for i in range(5):   # 5 sponges
        items.append({"accession": f"NC_9{i:05d}", "lineage": "Eukaryota; Metazoa; Porifera; Sp%d" % i, "n": i})
    for i in range(3):   # 3 cnidarians
        items.append({"accession": f"NC_8{i:05d}", "lineage": "Eukaryota; Metazoa; Cnidaria; Sp%d" % i, "n": i})
    score = lambda g: g["accession"]
    picked = balanced_pick(_build_trie(items, "Metazoa"), 12, score)
    assert len(picked) == 12
    counts = {}
    for g in picked:
        ph = g["lineage"].split("; ")[2]
        counts[ph] = counts.get(ph, 0) + 1
    assert counts["Porifera"] == 4 and counts["Cnidaria"] == 3 and counts["Chordata"] == 5
    again = balanced_pick(_build_trie(items, "Metazoa"), 12, score)
    assert [g["accession"] for g in again] == [g["accession"] for g in picked]

def test_support_parsing():
    assert _support("87.5/99") == {"alrt": 87.5, "ufboot": 99}
    assert _support("0.95") == {"support": 95.0}
    assert _support("100") == {"support": 100.0}
    assert _support(None) == {}

def test_root_tree_on_outgroup(tmp_path):
    nwk = tmp_path / "t.nwk"
    nwk.write_text("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
    taxa = [{"label": "A", "group": "Gammaproteobacteria"}, {"label": "B", "group": "Betaproteobacteria"},
            {"label": "C", "group": "Rickettsiales"}, {"label": "D", "group": "Jakobida"}, {"label": "E", "group": "Jakobida"}]
    t, how = root_tree(nwk, taxa, ["Gammaproteobacteria", "Betaproteobacteria"])
    assert "outgroup" in how
    kids = t.seed_node.child_nodes()
    sets = [sorted(n.taxon.label for n in k.leaf_iter()) for k in kids]
    assert ["A", "B"] in sets

def test_literature_trees_validate(tmp_path):
    root = Path(__file__).resolve().parents[1]
    trees = load_literature(root / "data" / "literature", tmp_path)
    assert {"rickettsiales_sister", "sister_to_alpha", "eukaryote_supergroups"} <= set(trees)
    for tid, t in trees.items():
        dt = dendropy.Tree.get(data=t["newick"], schema="newick", suppress_internal_node_taxa=True, preserve_underscores=True)
        tips = {n.taxon.label for n in dt.leaf_node_iter()}
        assert tips <= set(t["tips"]), tid
        for c in t.get("citations", []):
            pass
    refs = __import__("yaml").safe_load((root / "data" / "literature" / "references.yml").read_text())
    for t in trees.values():
        for c in t.get("citations", []):
            assert c in refs, c
