"""Offline end-to-end run of the pipeline on fixtures (quick profile)."""
import json
import pytest
from mitophy import extract_mito, fetch_bact, orthologs, sampling, align, trees, stats, provenance
from mitophy.site import build

pytestmark = pytest.mark.e2e

def test_pipeline_end_to_end(cfg, tmp_path):
    extract_mito.run(cfg)
    genomes = extract_mito.load_genomes(cfg)
    assert len(genomes) >= 20
    by = {g["organism"]: g for g in genomes}
    assert "Andalucia godoyi" in by and by["Andalucia godoyi"]["n_markers"] >= 30
    assert by["Homo sapiens"]["n_markers"] == 13 and by["Homo sapiens"]["transl_table"] == "2"
    fetch_bact.run(cfg)
    assert len(fetch_bact.load_assemblies(cfg)) == 8
    orthologs.run(cfg)
    hits = (cfg.workdir / "origin" / "hits.tsv").read_text().splitlines()
    assert len(hits) > 50
    sampling.run(cfg)
    sel = sampling.load_selected(cfg)
    assert 5 <= len(sel) <= 40
    align.run(cfg, "diversification"); align.run(cfg, "origin")
    for a in ("origin", "diversification"):
        m = json.loads((cfg.workdir / a / "align" / "manifest.json").read_text())
        assert m["n_taxa"] >= 5 and m["n_sites"] > 500
    trees.run(cfg, "diversification"); trees.run(cfg, "origin")
    for a in ("origin", "diversification"):
        tj = json.loads((cfg.resultsdir / a / "tree.json").read_text())
        assert tj["n_taxa"] >= 5 and "newick" in tj and tj["tips"]
    stats.run(cfg)
    assert (cfg.resultsdir / "diversification" / "charts" / "length_by_group.svg").exists()
    provenance.run(cfg)
    assert (cfg.resultsdir / "latest.json").exists()
    out = build.run(cfg, tmp_path / "_site")
    assert (out / "index.html").exists() and (out / "data" / "origin" / "tree.json").exists()
    html = (out / "origin.html").read_text()
    assert "Where do mitochondria branch" in html
