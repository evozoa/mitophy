import json, shutil
from pathlib import Path
import pytest
from mitophy.config import load_config

FIX = Path(__file__).parent / "fixtures"

@pytest.fixture(scope="session")
def cfg(tmp_path_factory):
    wd = tmp_path_factory.mktemp("work")
    c = load_config("quick", workdir=wd, resultsdir=wd / "results", threads=2)
    # point fetch-mito at the fixture file
    mito = c.sub("mito")
    shutil.copy(FIX / "mito_subset.gbff.gz", mito / "mito_subset.gbff.gz")
    (mito / "source.json").write_text(json.dumps({"files": {"mito_subset.gbff.gz": {"path": str(mito / "mito_subset.gbff.gz"), "last_modified": "fixture", "size": 0}}}))
    # restrict bacteria to the fixture proteomes and pre-place them (so fetch-bact needs no network)
    have = {p.stem for p in (FIX / "proteomes").glob("*.faa")}
    c.taxa_origin["bacteria"] = {g: [e for e in es if e.get("gcf") in have] for g, es in c.taxa_origin["bacteria"].items()}
    c.taxa_origin["bacteria"] = {g: es for g, es in c.taxa_origin["bacteria"].items() if es}
    pdir = c.sub("bact") / "proteomes"; pdir.mkdir(exist_ok=True)
    for p in (FIX / "proteomes").glob("*.faa"):
        shutil.copy(p, pdir / p.name)
    return c
