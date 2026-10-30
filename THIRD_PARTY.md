# Third-party resources

The resources below are not redistributed in this repository, and the licences in `LICENSE` and `LICENSE-DATA.md` do not apply to them.

- Oetjen et al. bone-marrow counts, GEO series GSE120221, and Sun et al. peripheral-blood counts, GEO series GSE128066. Sample accessions, URLs and hashes are in `reproduction/manifests/geo_original.json`.
- scIB human immune benchmark file: figshare, doi:10.6084/m9.figshare.12420968.v8, file 25717328. URL and hash are in `reproduction/manifests/artifacts.json`.
- CellTypist Immune_All_Low v2 model (Teichlab), https://celltypist.cog.sanger.ac.uk/models/Pan_Immune_CellTypist/v2/Immune_All_Low.pkl. Its SHA-256 and metadata are in `reproduction/results/gate/checkpoint.json`. Files that contain quantities computed from the model's coefficients, means and scales (`reproduction/results/main/tables/model_contributions.npz`, `reproduction/results/gate/model_only.npz` and the gene terms in `figures_src/data/`) remain subject to the model's terms.
- 10x Genomics PBMC3k matrix, used in the pilot: URL and hash in `reproduction/manifests/artifacts.json`.
- Inter typeface 4.000 by Rasmus Andersson, used for the figures, under the SIL Open Font License 1.1 (`figures_src/fonts/OFL.txt`).
