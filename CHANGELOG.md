# Changelog

Automatically generated entries for pipeline runs; hand-written entries for content changes.

## 2026-08-16 — pipeline run (profile `ci`)

- RefSeq mitochondrion release files: mitochondrion.1.genomic.gbff.gz (Thu, 09 Jul 2026 18:20:43 GMT); 17720 complete mitogenomes parsed
- Origin tree: IQ-TREE LG+F+G4, 103 taxa × 10199 sites, lnL -821164.8672, rooted on outgroup (Betaproteobacteria, Gammaproteobacteria, Bacteria)
- Diversification tree: FastTree LG GAMMA, 618 taxa × 4030 sites, lnL -2095922.65
- Sampled genomes: 618 (+0 / −0 vs previous run)
- Tools: python 3.12.13, biopython 1.88, pyhmmer 0.12.1, mafft 7.526, trimal 1.5.rev1, iqtree iqtree3 3.1.3, fasttree 2.2.0, datasets 18.35.0

## 2026-08-15 — pipeline run (profile `default`)

- RefSeq mitochondrion release files: mitochondrion.1.genomic.gbff.gz (Thu, 09 Jul 2026 18:20:43 GMT); 17720 complete mitogenomes parsed
- Origin tree: IQ-TREE LG+F+G4, 103 taxa × 10206 sites, lnL -821983.6695, rooted on outgroup (Betaproteobacteria, Gammaproteobacteria, Bacteria)
- Diversification tree: FastTree LG GAMMA, 618 taxa × 4017 sites, lnL -2071826.718
- Sampled genomes: 618 (+0 / −0 vs previous run)
- Tools: python 3.13.11, biopython 1.87, pyhmmer 0.12.2, mafft 7.505, trimal 1.5.rev1, iqtree iqtree3 3.1.1, fasttree 2.2.0, datasets 18.25.1
## 2026-08-15 — site created

- Initial release: literature hypotheses for the origin of mitochondria (4 trees), eukaryote supergroup tree with mitogenome traits, automated origin and diversification analyses, methods and provenance pages.

## 2026-08-16 — deep origin tree added

- `results/origin_deep`: LG+C20+F+G4 (PMSF guide LG+F+G4), full search, 1000 UFBoot + 1000 SH-aLRT, 103 taxa × 10,206 sites; run locally (`--profile deep`, 16 threads). Mitochondria monophyletic (100/100), sister to Rickettsiales (93.6/96).
