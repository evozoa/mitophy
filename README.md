# mitophy — the evolution of mitochondria

A static website (GitHub Pages) about the evolution of mitochondria that combines

1. **curated literature phylogenies** — hypotheses for the origin of mitochondria among Alphaproteobacteria and a
   supergroup-level view of mitochondrial genome diversification, with citations; and
2. **de novo analyses that update themselves** — a monthly GitHub Actions run downloads all RefSeq mitochondrial genomes
   and a curated set of alphaproteobacterial proteomes, extracts marker proteins, builds concatenated alignments,
   infers trees (IQ-TREE / FastTree), computes genome statistics and commits the results; the site is rebuilt from them.

Site: https://evozoa.github.io/mitophy/

## Layout

| path | purpose |
|---|---|
| `config/` | pipeline profiles, marker table (synonyms + bacterial homologues), taxon lists, sampling rules, site settings |
| `mitophy/` | Python package: pipeline stages (`fetch_mito`, `extract_mito`, `fetch_bact`, `orthologs`, `sampling`, `align`, `trees`, `stats`, `provenance`) and the site builder (`mitophy/site`) |
| `data/literature/` | hand-curated trees (Newick + JSON sidecar), hypotheses, references |
| `data/taxonomy/` | lineage → display-group/colour mapping |
| `static/` | CSS and JS (`mtree.js`, a dependency-free SVG tree viewer) |
| `results/` | committed pipeline outputs (trees, alignments, tables, charts, run records) |
| `work/` | intermediates (git-ignored) |
| `tests/` | unit tests and an offline end-to-end fixture run |
| `.github/workflows/` | `pipeline.yml` (monthly run + deploy), `pages.yml` (deploy on content changes), `test.yml` |

## Quick start

```bash
micromamba create -y -f environment.yml -n mitophy && micromamba activate mitophy   # or: pip install -e . (+ mafft, trimal, iqtree, fasttree, ncbi-datasets-cli on PATH)
mitophy all --profile quick --workdir work-quick     # smoke run (~10 min incl. RefSeq download)
mitophy all                                          # full local run (origin tree with UFBoot: ~1 h on 16 cores)
mitophy all --profile deep --threads 20              # site-heterogeneous origin tree, published to results/origin_deep
mitophy site && python -m http.server -d _site 8000  # build and preview the site
pytest -q && pytest -q -m e2e                        # tests
python scripts/refresh_bacterial_accessions.py       # re-resolve bacterial taxa -> RefSeq assemblies
python scripts/check_dois.py                         # verify all DOIs resolve
```

Stages: `fetch-mito → extract-mito → fetch-bact → orthologs → sample → align → trees → stats → provenance → site`.
Every stage writes a `manifest.json` with an input hash and is skipped when nothing changed (`--force` overrides).

## Licence

Code: MIT. Site content (text, curated trees): CC-BY 4.0. Data from NCBI RefSeq is in the public domain.
