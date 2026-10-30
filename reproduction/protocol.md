# Query composition and zero-gene filtering in released CellTypist annotation

Protocol version: 0.4 draft, 23 October 2026. Status: prospective main-study design, pending feasibility gate and external timestamp. Pilot analyses and subsequent amendments will be disclosed. The exploratory PBMC3k gate pilot has now been run and is disclosed in results/gate/gate_report.md. Main-study annotation endpoints have not been computed.

## 1. Objective and scope

Determine whether removing genes recorded as zero throughout a query changes a released CellTypist model's predictions, and whether independently filtering different query compositions changes the prediction for an identical cell. Test whether these changes follow the model's linear scores and whether restoring the recorded feature schema recovers the reference output.

Two questions are distinct: retention versus deletion isolates the representation change; pooled versus per-sample filtering tests composition dependence. Differences do not establish which label is biologically correct. Confidence is the native model score, not an assumed calibrated probability of correctness.

The proposed contribution is a reproducible empirical characterization of this specific preprocessing interaction and its practical consequences. Linear-model algebra, missing-feature handling and mean centring are established methods. The study does not claim a new classifier, a new calibration method or priority established by absence of search results.

## 2. Model and mechanistic prediction

Use released `Immune_All_Low` v2, best-match mode, with majority voting disabled. Record the exact official source, file checksum, model metadata, feature/class order, scaler parameters and software environment. If the named version is unavailable, do not silently use a different checkpoint.

For an observed-zero model gene g, define a_kg = w_kg * mu_g / sigma_g, where k indexes classes and parameters come from the saved centred scaler and classifier. For an omitted zero-gene set Z:

    delta_k = sum over g in Z of a_kg
    z_B,ik = z_A,ik + delta_k
    predicted_label_B,i = argmax over k of z_B,ik
    confidence_B,i = sigmoid(max over k of z_B,ik)

Verify the actual implementation's feature alignment, intercept, class conventions, scaling and clipping. For nonnegative expression/training means, standardized zeros are nonpositive and the documented upper clipping does not alter them. For an uncentred model this zero-column contribution is zero; verify the expected null and stop the centred-checkpoint route. Do not change the checkpoint to manufacture an effect.

For arbitrary removed features in secondary arm D, use delta_ik = -sum over removed g of w_kg * scaled_clipped(x_ig). This shift can vary by cell.

Before expression data, save every a_kg and two gene diagnostics: R_g = max_k(a_kg) - min_k(a_kg), and M_g = max_k(abs(a_kg)). Rank by R_g, breaking ties by gene name. The spread of the summed shift is bounded by sum(R_g) over omitted genes. A unique winner whose original margin over its runner-up exceeds that spread is stable in exact arithmetic. The converse does not hold. Evaluate all class scores to predict actual flips. Common class shifts can alter sigmoid confidence without changing labels.

These diagnostics are conditional mathematical quantities, not gene importance, deletion probabilities or evidence of a practically attainable worst case. They must not influence masks, sample selection, QC or partition seeds.

## 3. Data eligibility and predetermined fallback

Pilot: original PBMC3k 10x count matrix, barcodes and feature table, retaining original gene IDs and symbols. Do not use the convenience H5AD as the sole identifier source: its documented preparation makes symbols unique. It is an implementation resource, not independent main-study evidence. Artificial partitions are not natural samples.

Main input: original GEO count matrices for the three Oetjen and four Sun natural samples, restricted to the same 18,410 cells selected in the scIB human immune benchmark. Recover each cell by its original barcode and sample suffix; require exact equality of GEO and scIB counts for every one of the 12,303 shared genes in every retained cell before eligibility. Retain all original recorded feature rows for normalization, then apply the fixed identifier exclusions. GEO accessions, file hashes, cell mappings and prepared-input hashes are saved in manifests. Do not apply new biological QC. The reduced scIB count matrix is retained only as a labelled secondary input sensitivity, not independent evidence. Its pooled all-zero filtering and cross-study intersection are explicitly documented in amendments.md. This provenance amendment follows PBMC3k pilot inspection but precedes main annotation outcomes.

**The first data check is the sample table.** For every recorded batch, establish source study, acquisition sample, donor if known, tissue/technology, number of cells, expression representation and measured feature support. A batch column alone is not a sample/donor definition. Unknown donor identities remain unknown. Technical libraries may be computational query units when their identity is documented, but they do not become independent biological replicates.

For the practical primary endpoint, each eligible study must have at least two documented natural acquisition samples using a compatible measured gene schema. The parent pool is all eligible samples within that study. A study containing only one eligible sample contributes to the mechanism endpoint only.

Scope is selected from provenance, never from effect size:

| Eligible studies with at least two compatible natural samples | Prespecified scope |
|---|---|
| Two or more | Limited cross-study technical note; no population-wide guarantee |
| One | Accept a single-study composition case study; remove cross-study practical claims |
| Zero | Mechanism/reproducer report only, not a preprint under this plan; random partitions cannot rescue the natural-sample claim |

The GEO repair retrieves original counts for the already selected samples, rather than selecting a new cohort. No automatic backup-dataset search is included in the study. Recovery of documented metadata for the same candidate is allowed. A new dataset requires a dated amendment and new pre-analysis freeze before its predictions; adding data after seeing results is disclosed as an extension.

Whole studies and samples are excluded only for documented input/provenance failures, not because they have weak effects. Record exclusions and counts. Two eligible studies are a scope minimum, not a power calculation. Restricted gene sets are not assumed to underestimate effects: signed contributions and margins can change sensitivity in either direction.

## 4. Expression preparation and reference definition

Use documented nonnegative counts or documented compatible log1p CP10k values. A layer named `counts`, integer-like values or `.raw` is not sufficient proof. Establish provenance per study; mixtures of counts, TPM or other transforms cannot be treated as uniform counts. Do not reconstruct counts by rounding, invert undocumented normalization or use integrated/standardized expression as counts.

Preserve the supplied cell QC. Exclude empty-library cells from count normalization; report them. Resolve non-finite/negative inputs before a sample is eligible.

**Fixed duplicate policy.** Preserve original stable feature IDs, symbols and row identities. Identify symbols assigned to multiple feature rows within any source feature table. For each study parent, take the union of those ambiguous symbols across its eligible sources; exclude every row carrying any such symbol from scoring in every sample, pool and arm, even if another sample has only one copy. The same gene occurring once in each sample is not a duplicate. Exclude individually identified missing/unresolved symbols as well. Define and save this exclusion from feature metadata before query partitioning, without using expression or predictions. Do not select the first copy, aggregate distinct features, silently suffix names, or strip suffixes from apparently renamed names. Genuine hyphenated symbols must not be altered.

If an artifact already has renamed features, recover the exact original identity from retained stable IDs or traceable original feature tables. If this correspondence cannot be established, that source is ineligible; unique-looking names alone do not clear the check. Log every excluded row, original ID/symbol, reason and model membership, plus feature overlap and removed expression mass. Assert that all remaining scoring symbols are unique. These exclusions are fixed preprocessing, not B/D intervention genes.

Normalize verified counts once to 10,000 over all recorded available features, including those subsequently excluded for identifier ambiguity, then log1p. Physically remove the fixed excluded rows from the scoring matrix only after normalization. If this is a reduced-feature count matrix, explicitly call this normalization over available features, not reconstruction of the original full transcriptome. Use already compatible normalized values without a second normalization. Do not renormalize after fixed exclusions or in any filtered arm. Keep the original feature table and per-cell denominator/normalization provenance.

Identify actual measured zeros separately from artificial zeros introduced by concatenation. Establish one compatible schema within each study. If a verified common schema is required, define it before predictions and describe the resulting restriction. Do not fill unmeasured or unknown features with zeros.

**Arm A is the reference because it preserves the recorded measurements within the fixed unambiguous scoring schema.** Under the checked centred model, an explicit zero receives its zero-expression contribution; omitting it removes that contribution, algebraically equivalent to setting that feature to its training mean for scoring. This is a controlled computational baseline, not proof that A is biologically correct or that its complete input exactly matches the model's training data. All results are conditional on the common identifier exclusions; their effect is not assumed negligible.

Provide each arm through an explicit matrix. Record the matrix and features actually consumed by CellTypist. Prevent `.raw` or another layer from silently restoring a different schema. Excluding expressed features can make the inverse-log row sum less than 10,000. Record and reconcile this deficit; an explained sum-check warning is not an input rejection and does not justify renormalization. Inspect the pinned implementation's actual validation behavior. Unexpected rejection or substitution blocks the affected comparison.

## 5. Arms, query contexts and computation

| Arm | Intervention |
|---|---|
| A | Preserve the fixed unambiguous parent scoring schema |
| B | Delete only columns exactly zero throughout that semantic query |
| C | Restore exactly the recorded zero columns deleted in B |
| D, secondary | Delete genes detected in fewer than three query cells (`min_cells=3`) |

Detection means positive expression on verified nonnegative data. D removes some nonzero information; record its removed entries and expression and report it separately. It is a fixed-normalization filtering experiment, not a universal reconstruction of workflows that filter before normalization. Do not recover D's removed nonzero entries by zero filling.

For each main study, evaluate A/B/C on each natural sample and on the within-study pool. Compare identical cell IDs across contexts. D repeats these contexts as a secondary analysis. Cross-study pooling is omitted from the one-day main plan to avoid introducing another feature-support comparison.

Stress test: random partitions within each eligible study, at sizes 32, 128 and 512, for seeds 0, 1, 2, 3 and 4. Use a seeded permutation of sorted globally unique cell IDs; retain and report final smaller chunks. Save actual membership files so reproduction does not depend solely on RNG implementations. Evaluate A/B/C; D is restricted to natural sample/pool contexts. These results measure computational sensitivity, not extra biological replication.

Distinguish a semantic query from a computation chunk. Compute a query's mask using all its cells, then keep that mask fixed while scoring in memory-bounded chunks (default 512 cells). Changing computation chunk size must not change the mask. A single CellTypist API call is not required to contain an entire semantic query. Start each arm/query from a pristine model instance or verified independent copy.

## 6. Endpoints and summaries

Two primary endpoints are reported separately:

1. **Mechanism:** per natural sample u, the fraction of identical cells with different native labels in A_u versus B_u.
2. **Practical composition:** per sample u, the fraction with different labels in B_u versus independently filtered B_pool from its within-study parent pool.

Report every sample's originally eligible count, valid paired count, changed count and failures. Fractions use valid paired predictions; incomplete comparisons are prominently flagged and cannot support a claim of complete primary results. Failed predictions, zero-feature inputs, NaNs and unmatched IDs are not unchanged labels. Do not exclude numerical near ties from the principal denominator.

Within-study summaries are total changed cells divided by total valid paired cells. Across-study summaries are the equal-weight mean of those study fractions; show eligible study sets separately for each endpoint. Do not use a multi-study summary for a one-study case. Report sample-specific results alongside aggregates.

Secondary measurements: per-query class shift vectors, per-A-label change counts/rates for labels represented by at least 50 cells within a natural sample, and counts of cells certified stable by their top-two margin exceeding the spread of the exact shift. For the composition contrast, certify stability using pooled B scores and the sample-minus-pool shift. These additions were selected after the pilot but before main predictions. No within-lineage or neighbouring-subtype claim is planned without a separately frozen lineage map. Additional secondary measurements: all-class score residuals; predicted versus observed destinations and flip indicators; mean/median/quantiles of confidence change; threshold crossings in both directions at 0.50, 0.90 and 0.95 (acceptance means >= threshold); different labels with confidence >=0.90 in both compared arms; predicted composition changes; numbers and contributions of removed model genes; all stress-test settings; D results.

Report confidence and composition changes for both primary contrasts, not just cells that flip. Fixed-threshold changes do not prove calibration or correctness. No biological accuracy benchmark, calibration fitting, new classifier training, bootstrap confidence intervals or confirmatory p-values is planned for this one-day version. The outputs are exact descriptive results on the included finite datasets. These limits avoid treating cells, technical samples or repeated partitions as independent population evidence.

## 7. Required validation

Predict B scores, labels, confidence and threshold decisions from A plus the analytic shift before running the independent B API calculation. Use all classes and the checkpoint's class order. Use the analogous per-cell reconstruction for D. A common shift, cancellation and a non-runner-up overtaking the winner must be included in synthetic checks.

Default validation: score absolute tolerance 1e-6 and relative tolerance 1e-6; floating-point diagnostic calculations use float64. For each cell define epsilon from the absolute tolerance plus relative tolerance times the maximum absolute compared score; a top-two margin <=2*epsilon is flagged as numerically ambiguous. Save exact API labels and all mismatches. Tolerances may be revised only through documented synthetic/pilot checks before main-data predictions, with reasons and original results retained. Never relax them to excuse an unexplained main-study failure.

Required controls: celltypist.annotate(mode="best match", majority_voting=False) versus the direct Classifier path on the PBMC3k pilot; gene-order permutation; fixed-schema pooled versus partitioned scoring; deletion of all-zero non-model genes; C recovery of A; direct score reconstruction versus API output; repeated runs with a freshly loaded checkpoint and reversed arm order; bounded dense/sparse equivalence. For label invariance require agreement outside recorded near ties, and preserve all tie-related changes in raw results. Check stable ID joins and equality of retained expression entries.

Identifier tests must include two rows sharing one symbol with one zero and one expressed: both are excluded in every arm, with the original normalization denominator preserved. Also test that row permutation cannot select a different copy, that study-wide exclusions remain identical across sample/pool contexts, and that genuine hyphenated names are not stripped.

Unexpected primary residuals block interpretation until explained. Document explanations and rerun affected controls. Model-only reconstructed results cannot replace a failed primary API validation. D failure from the input contract is reported separately; do not silently renormalize, bypass validation or let `.raw` replace its matrix. Secondary D failure does not invalidate successful A/B/C results.

## 8. Prospective private snapshot and amendments

The PBMC3k pilot preceded main inference and is disclosed. Version 0.4 replaces the earlier public-before-main rule after the author explicitly requested that work remain private while the project be completed. This privacy-driven amendment is recorded before main annotation outcomes. Preserve the earlier v0.3 bundle and its hashes as history.

Before main inference, save an allowlisted snapshot of the protocol, amendments, fixed and resolved configuration, source, tests, dependency lock and relevant input/query manifests. Hash every file and the snapshot archive, record the actual UTC creation time and the author's private-completion instruction, and verify all source/input hashes before the run. The private snapshot supplies internal integrity and execution provenance; its timestamp is not independently verifiable public preregistration. No external timestamp has been obtained. Do not represent local dates, private commits or later publication as prospective public evidence.

Main inference is permitted under this explicit private-snapshot amendment when technical inputs and controls pass. All outcomes remain private until the author reviews and authorizes release. At release, publish the preserved snapshot, pilot disclosure, execution logs suitable for sharing, input retrieval hashes and all amendments. Any later implementation repair records whether outcomes had been observed, retains previous files, and reruns affected checks. Changes to scientific endpoints after outcomes are exploratory or require a separately disclosed extension.

## 9. Interpretation, deliverables and release criteria

Prepare a complete report of all valid planned results irrespective of effect size, with public release controlled by the author. A failed primary control blocks scientific release until resolved. A completed, validated two-study experiment may support a technical preprint, including a null result; acceptance or novelty is not guaranteed. All limitations must be explicit and claims must match the eligible scope. For headline wording only, define a material effect as at least 1.0% label disagreement among valid paired cells in at least one eligible natural sample on either primary endpoint. This threshold was explicitly fixed after the PBMC3k pilot, before main outcomes; do not claim it was prespecified before the pilot. It is not a significance test, clinical threshold, selection rule or condition for reporting valid results. A mechanism effect does not establish natural-sample composition dependence. Artificial-chunk findings alone support only the computational sensitivity claim. A null does not justify universal stability claims. All planned results, failures and amendments remain in the report irrespective of publication value.

Deliver three figures: natural-sample primary outcomes; analytic versus observed changes with residuals and label checks; random-partition sensitivity with fixed-schema controls. Deliver two main tables: sample/input provenance and primary outcomes. Supplement D, threshold crossings, feature contributions, validation and deviations. Save complete per-cell predictions/scores and machine-readable summaries so figures and manuscript numbers can be regenerated without rerunning inference.

Use neutral language about intended behavior; seek maintainer clarification with a minimal reproducer if authorized. Maintainer response is not required for numerical validation and silence is not endorsement. Independently review all manuscript claims and citations. Acknowledgment, if accurate at release: “AI tools assisted with code development and manuscript drafting; the authors reviewed and take responsibility for the methods and results.” Record the actual tools and roles used. Public posting and screening time are not guaranteed by an analysis deadline.

## Sources

- [CellTypist inference implementation](https://celltypist.readthedocs.io/en/stable/_modules/celltypist/classifier.html) and [FAQ](https://www.celltypist.org/faq).
- [scTab](https://www.nature.com/articles/s41467-024-51059-5), including its centring choices; [CellTypist issue 29](https://github.com/Teichlab/celltypist/issues/29).
- [scIB benchmark](https://www.nature.com/articles/s41592-021-01336-8).
- [Open Problems immune dataset listing](https://openproblems.netlify.app/datasets/openproblems_v1/immune_cells/). Listed slots are metadata, not per-study input validation.
- [PBMC3k access documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html).
- [Scanpy dataset source](https://github.com/scverse/scanpy/blob/main/src/scanpy/datasets/_datasets.py) documents preparation of the convenience PBMC3k artifact, including unique-name conversion. Record the installed version and original 10x artifact at execution.
- [scIB immune preparation notebooks at inspected revision](https://github.com/theislab/scib-reproducibility/tree/3afbffd3674726e5146797be21cf6bd7470a2c5f/notebooks/data_preprocessing/immune_cells). Verify historical concatenation behavior and the selected artifact's full derivation; the repository alone does not establish measured-zero provenance for every derivative.
- [CPA immune-data example](https://cpa-tools.readthedocs.io/en/stable/tutorials/Batch_correction_in_expression_space.html), which notes missing counts for Villani in its artifact; verify applicability to the selected file.
- [OSF registration guidance](https://help.osf.io/article/330-welcome-to-registrations), distinguishing a timestamped read-only study plan from ordinary project preparation.
