"""Stage 3: download RefSeq proteomes for the bacterial taxa in config/taxa_origin.yml via NCBI `datasets`."""
from __future__ import annotations
import csv, json, logging, shutil, zipfile
from pathlib import Path
from .config import Config
from .cache import hash_inputs, up_to_date, write_manifest
from .tools import run as sh
from .taxonomy import GroupMapper

log = logging.getLogger("mitophy")

def bacterial_entries(cfg: Config) -> list[dict]:
    out = []
    for group, entries in cfg.taxa_origin["bacteria"].items():
        for e in entries:
            if e.get("gcf"):
                out.append({**e, "config_group": group})
    maxb = cfg.params["origin"].get("max_bact")
    if maxb:
        # keep group diversity in quick mode: round-robin over groups
        by_group: dict[str, list] = {}
        for e in out:
            by_group.setdefault(e["config_group"], []).append(e)
        picked, i = [], 0
        while len(picked) < maxb and any(by_group.values()):
            for g in list(by_group):
                if by_group[g] and len(picked) < maxb:
                    picked.append(by_group[g].pop(0))
        out = picked
    return out

def run(cfg: Config) -> Path:
    outdir = cfg.sub("bact")
    entries = bacterial_entries(cfg)
    accs = sorted(e["gcf"] for e in entries)
    ih = hash_inputs(objs=[accs])
    table = outdir / "assemblies.tsv"
    if up_to_date(outdir, ih, [table] + [outdir / "proteomes" / f"{a}.faa" for a in accs], cfg.force):
        log.info("fetch-bact up to date (%d proteomes)", len(accs))
        return table
    pdir = outdir / "proteomes"; pdir.mkdir(exist_ok=True)
    missing = [a for a in accs if not (pdir / f"{a}.faa").exists() or cfg.force]
    if missing:
        acc_file = outdir / "accessions.txt"; acc_file.write_text("\n".join(missing) + "\n")
        zip_path = outdir / "ncbi_dataset.zip"
        dbin = cfg.tools.get("datasets", {}).get("bin", "datasets")
        sh([dbin, "download", "genome", "accession", "--inputfile", acc_file, "--include", "protein",
            "--filename", zip_path, "--no-progressbar"], log_file=outdir / "datasets.log")
        with zipfile.ZipFile(zip_path) as z:
            for name in z.namelist():
                parts = name.split("/")
                if len(parts) >= 4 and parts[1] == "data" and name.endswith("protein.faa"):
                    acc = parts[2]
                    with z.open(name) as src, open(pdir / f"{acc}.faa", "wb") as dst:
                        shutil.copyfileobj(src, dst)
        zip_path.unlink()
        still = [a for a in missing if not (pdir / f"{a}.faa").exists()]
        if still:
            log.warning("no protein.faa retrieved for: %s", ", ".join(still))
    mapper = GroupMapper(cfg.supergroups)
    # lineage: use datasets summary if available in cache, else derive from config group
    with open(table, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "name", "gcf_name", "config_group", "group", "supergroup", "level", "n_proteins"])
        for e in entries:
            faa = pdir / f"{e['gcf']}.faa"
            n = sum(1 for l in open(faa) if l.startswith(">")) if faa.exists() else 0
            grp = mapper.group_of([e["config_group"], "Alphaproteobacteria" if e["config_group"] != "Outgroups" else "Bacteria"], e["name"])
            if e["config_group"] == "Outgroups":
                grp = _outgroup_group(mapper, e["name"])
            w.writerow([e["gcf"], e["name"], e.get("gcf_name", ""), e["config_group"], grp["id"], grp.get("supergroup", ""), e.get("level", ""), n])
    write_manifest(outdir, ih, n_proteomes=len(accs))
    return table

_BETA = {"Neisseria", "Burkholderia", "Nitrosomonas", "Ralstonia", "Cupriavidus", "Bordetella"}
def _outgroup_group(mapper: GroupMapper, name: str) -> dict:
    genus = name.split()[0]
    return mapper.group_of(["Bacteria", "Betaproteobacteria" if genus in _BETA else "Gammaproteobacteria"], name)

def load_assemblies(cfg: Config) -> list[dict]:
    p = cfg.workdir / "bact" / "assemblies.tsv"
    if not p.exists():
        raise SystemExit("run `mitophy fetch-bact` first")
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
