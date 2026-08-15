#!/usr/bin/env python
"""Resolve bacterial taxa in config/taxa_origin.yml to RefSeq assembly accessions (GCF_) using `datasets summary`.
Preference: reference genome > complete genome > chromosome > any RefSeq assembly. Writes `gcf:` and `gcf_name:` back
into the YAML (idempotent; entries with `gcf` are skipped unless --force)."""
import argparse, json, subprocess, sys
from pathlib import Path
import yaml

LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}

def summary(taxon: str, extra: list[str]) -> list[dict]:
    cmd = ["datasets", "summary", "genome", "taxon", taxon, "--assembly-source", "RefSeq", "--as-json-lines", *extra]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = []
    for line in p.stdout.splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    return out

def pick(recs: list[dict]) -> dict | None:
    def key(r):
        ai = r.get("assembly_info", {})
        cat = ai.get("refseq_category", "")
        return (0 if "reference" in cat else 1 if "representative" in cat else 2,
                LEVEL_RANK.get(ai.get("assembly_level"), 9),
                -int(r.get("assembly_stats", {}).get("total_sequence_length", 0)))
    return min(recs, key=key) if recs else None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/taxa_origin.yml")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    path = Path(a.config)
    cfg = yaml.safe_load(path.read_text())
    n_ok = n_fail = 0
    for group, entries in cfg["bacteria"].items():
        for e in entries:
            if e.get("gcf") and not a.force:
                continue
            taxon = e.get("taxon", e["name"])
            recs = summary(taxon, ["--reference"]) or summary(taxon, ["--assembly-level", "complete,chromosome"]) or summary(taxon, [])
            best = pick(recs)
            if best:
                e["gcf"] = best["accession"]
                e["gcf_name"] = best["organism"]["organism_name"]
                e["level"] = best.get("assembly_info", {}).get("assembly_level")
                n_ok += 1
                print(f"[ok]   {group:28s} {e['name']:40s} {e['gcf']}  {e['gcf_name']} ({e['level']})")
            else:
                n_fail += 1
                print(f"[FAIL] {group:28s} {e['name']:40s} no RefSeq assembly for taxon '{taxon}'", file=sys.stderr)
    header = "".join(l for l in path.read_text().splitlines(keepends=True) if l.startswith("#"))
    path.write_text(header + yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True, width=200))
    print(f"resolved {n_ok}, failed {n_fail}")

if __name__ == "__main__":
    main()
