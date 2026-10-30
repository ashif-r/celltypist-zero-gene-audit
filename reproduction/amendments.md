# Prospective amendments and disclosure

## Version 0.3 — 23 October 2026, before main annotation outcomes

The PBMC3k pilot was already known: 86/2,700 labels differed in A versus B and the analytic formula matched. Main annotation endpoints had not been computed. Only input metadata, gene tables, raw counts, masks and barcode/count reconciliation were inspected for the main samples. The feedback reports an independent pilot check, but that checker was not supplied as an executable file in this workspace; we do not claim to have run it.

1. **Input provenance repair.** The scIB Oetjen and Sun merge notebooks call `filter_genes(min_cells=1)` before cross-study gene intersection. The earlier gate documented the intersection but missed its consequence: study-pool B and A become identical on this reduced schema, so the old composition and mechanism comparisons are redundant. Input-only checks found zero all-zero genes in both pooled scIB studies. Replace main expression with original GEO counts, retaining the same seven samples, 18,410 cells and metadata. Verify all 226,498,230 shared gene/cell count entries exactly; record all source files, barcodes and hashes. Normalize over all original recorded features before excluding duplicate symbols. Keep the scIB input only as a labelled secondary sensitivity. This change is motivated by provenance and precedes main outcomes; it is not an effect-driven dataset search.

2. **Secondary diagnostics added after pilot inspection.** Report class-specific query shift vectors, per-A-label changes for labels with >=50 cells in a natural sample, and margin-versus-shift-spread stability certificates. Composition certificates use pooled-B reference scores and the sample-minus-pool shift. These are descriptive and pilot-informed. No unsupported neighbouring-subtype/lineage claim is planned.

3. **Materiality and reporting.** Fix >=1.0% label disagreement among valid paired cells in any eligible natural sample on either primary endpoint as a headline-wording rule. No reliable evidence was found that the authoritative v0.2 config fixed this threshold before the pilot. Disclose its post-pilot selection. Report all valid planned results regardless of effect size. Failed primary controls block release; no promised bioRxiv acceptance follows from this policy.

4. **Implementation checks.** Add public annotate-versus-Classifier equality, prospective query/mask manifests, substantive mask/partition/failure tests, and integration tests using the actual native API. Preserve original gate results as historical evidence.

5. **Privacy and external freeze.** The current private repository remains private. A clean allowlisted release folder will be prepared with protocol, amendments, configuration, code, tests, input provenance/checksums and relevant manifests. It excludes private/, expression downloads, model files and logs. No main annotation inference is authorized by an absent freeze record. The public-readable rule of v0.2 remains in force until the author chooses and authorizes either that route or a specifically disclosed alternative prospective timestamp amendment. A local/private commit alone is not an external timestamp.

## Version 0.4 — 23 October 2026, before main annotation outcomes

The author explicitly directed: “keep it private - i will go through it first. but, complete the project.” To honor this direction, replace the public-readable/external-before-main rule with a private prospective source/input snapshot and recorded hashes, verified before main inference. This amendment follows the known PBMC3k pilot, repaired input counts and masks, but precedes main annotations. No main endpoints were computed under the old rule.

This is an author-directed privacy amendment, not an external registration. Local UTC dates and private GitHub commits are not independent timing evidence. Disclose that no independent public preregistration was obtained. Preserve v0.3 and the exact private v0.4 bundle and release them with the eventual paper only after author approval. The numerical design, repaired cohort, endpoints, tests, materiality wording and requirement to report all valid results are unchanged.

## Implementation revision 0.4.1 - 2026-10-27T04:43:18.034108+00:00

The first main run was deliberately interrupted for performance after 79 completed queries. Those annotation outputs were available; this revision is not prospective relative to them. Preserve the entire partial run and original v0.4 snapshot. The scientific protocol, cohort, normalization, masks, query memberships, endpoints, tolerances and reporting rule remain unchanged.

Before native calls, validate the full matrix and omit only non-model features that the native API itself ignores; do not renormalize. Replace object-array membership operations with ordered set lookups. Use four isolated processes, fresh model copies for each native call, bounded query submission and deterministic saved output order. Check exact equivalence against the preserved native outputs and an eight-query isolated-worker control before the full restart. The full final run will be independently audited and repeated from clean output directories. This is an implementation repair after partial outcomes, not a new preregistration or an external timestamp.
