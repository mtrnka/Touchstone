# Touchstone stabilization roadmap

This roadmap records the agreed direction for making the current package a
stable successor to the original script-based Touchstone workflow. The public
data model will remain a flat data frame. More complicated internal structures
may be used where helpful, but package functions should return canonical flat
tables with enough provenance to understand how each result was produced.

The core binary-crosslink MS2 workflow is now substantially stabilized and has
been exercised on TRiC, translocon HCD/EThcD, large E. coli, and Astral data.
The full Astral peptide-pair analysis and report have been delivered, and the
working ribosome example and GitHub README are complete. The remaining release
priority is therefore a clean-install test, complete package check, and stable
public MS2 tag. PPI-specific scoring, MS3 reconstruction, and deeper Prospector
integration remain important but are not blockers for that release.

Unless an analysis explicitly states otherwise, diagnostics, thresholds, FDR
summaries, and reported counts must keep intra-protein and inter-protein results
separate at every summarization level.

## 1. Stabilize the binary-crosslink workflow

- [x] Repair package loading, dependency declarations, and generated
  documentation.
- [x] Isolate Protein Prospector Search Compare column normalization and test
  common report-column selections.
- [x] Add small, hand-checkable tests for decoy/FDR calculations, pair
  construction and summarization, and product-ion evidence.
- [x] Make model fitting reproducible: control randomness and prevent spectra
  or equivalent crosslinks from leaking across training and validation groups.
- [x] Define conservative feature profiles for approximately small (up to
  10--20 proteins), medium (20--200), and large (more than 200) systems. Treat
  these boundaries as starting defaults, not biological laws, and report the
  chosen profile. Automatic selection counts both accessions participating in
  plausible high-Score.Diff intra- and inter-protein target CSMs without
  requiring repeated observations of the same URP, and retains the original
  dominant-protein safeguard for small systems with background matches.
- [x] Add explicit training and result objects, model-selection diagnostics,
  score-correlation plots, and reproducible decoy downsampling while retaining
  direct access to all fitted candidates.
- [x] Implement and evaluate training-only `Score.Diff` prefiltering separately
  from feature-complexity selection. The prefilter limits model fitting without
  removing lower-Score.Diff CSMs from subsequent scoring. Product-ion evidence
  remains a reporting-polish option rather than a general training prefilter.
  A fixed-threshold override supports reproducible validation and comparison
  with historical analyses. The legacy low-FDR interprotein summary
  (`interInt`) remains internal to automatic prefilter selection and is not
  exposed in the hyperparameter candidate audit.
- [x] Make a linear SVM the default model. When explicitly requested, radial
  candidates receive a separate recommendation while the overall default
  remains linear. The validation datasets did not establish a reproducible
  benefit for radial models, and several radial results were visibly or
  structurally less credible, so radial support remains exploratory.
- [x] Complete the candidate audit. The candidate table records filters,
  features, parameters, validation FDR/yield, relative recovery within each
  kernel family, concise intra- and interprotein rank correlations over the full
  range and high-Score.Diff tail, eligibility, recommendation status, and the
  selection or rejection reason. Poorly fitted candidates are excluded before
  the near-best recovery comparison, and selection falls back to intraprotein
  recovery when credible models recover no interprotein links.
- [x] Measure repeat-split stability and make a three-fit score ensemble the
  default for the selected recommendation. Prefilter selection and the
  hyperparameter grid run once; only the selected linear model (and an
  explicitly requested radial recommendation) is refit with consecutive seeds.
  CSM scores are averaged before polishing, summarization, and thresholding.
  The ensemble closely reproduced the two-of-three consensus on translocon HCD,
  E. coli URPs/PPIs, and Astral DSBSO peptide pairs. Set
  `ensembleRepeats = 1` for the former single-fit behavior.
- [x] Stabilize threshold estimation. Retain modeled interprotein thresholds to
  handle sparse, jagged empirical FDR curves, with documented fallbacks when a
  fit fails; do not interpret the raw count-based `calculateFDR()` value as
  closer to ground truth when decoys are sparse or impose an automatic rule
  that the raw estimate may never exceed the requested target. Keep raw counts
  available as diagnostics without overriding the modeled estimate.
- [x] Separate statistical classification from evidence polishing. Provide
  named, transparent polishing policies (for example, minimum product ions or
  backbone-ladder coverage), record the applied rules and their before/after
  CSM counts, and recalculate summaries on the polished result. Three distinct
  product-ion cleavages per peptide are the current general default when the
  required ion annotations are available; more stringent ladder policies remain
  explicit, dataset-specific reporting choices.
- [x] Add the classified reporting path: apply the reporting-level thresholds
  to scored CSMs, rerun `calculatePairs()` on the surviving CSMs, and only then
  summarize to URP, peptide-pair, protein-pair, or module-pair. This ensures
  `numCSM` and `numURP` describe threshold-passing evidence; the training-only
  weighted support features are omitted from classified reporting tables.
- [x] Run the MS2 testing protocol across representative datasets spanning
  small and large systems, HCD and EThcD fragmentation, DSSO and DSBSO
  crosslinkers, and Orbitrap/Astral acquisition. The local harness and saved
  audits cover TRiC, translocon, E. coli with a human entrapment database, and
  Astral DSBSO. Future datasets may extend this suite, but additional examples
  are not required before the public MS2 release.
- [x] Validate the current small-system workflow against the collaborator-
  delivered TRiC sample-1 result. The historical thresholding and ladder-
  polishing order reproduced the same total of 554 target URPs, sharing 277 of
  278 inter-protein and 272 of 276 intra-protein URPs. Retain the historical
  sequential-ladder requirement as an explicit conservative reporting policy;
  do not use the high open/closed-structure violation rates as a truth set for
  further tuning.
- [x] Validate HCD and EThcD translocon results sufficiently for stabilization.
  Peptide-2 fragmentation and distinct-ion inspection explained much of the
  questionable low-score EThcD recovery, and the minimum-three-ion policy
  improved the structural comparison. Ribosomal structural distances remain an
  imperfect truth set because conformational heterogeneity, polysomes, and
  aggregation can produce real over-length links. Preserve the analysis and
  conclusions without further structure-mapping work at this stage.
- [x] Re-test the training-only Score.Diff prefilter on the large E. coli data
  used to develop the original procedure. Fixed linear models at Score.Diff
  thresholds 0, 5, 10, 15, and 20 were scored against the complete input, with
  the historical minimum-four-product-ions-per-peptide rule applied before
  independent 2% CSM and protein-pair thresholds. Thresholds 5--20 produced a
  stable plateau, and the automatic value of 10 closely reproduced the
  historical CSM and PPI recovery while retaining lower-Score.Diff CSMs during
  full-dataset scoring.
- [x] Replace the self-inclusive, hard-thresholded `wtCSM` and `wtURP` features
  with self-excluded, bounded `CSMsupport` and `URPsupport`. Strengthening
  `URPsupport` did not improve E. coli PPI recovery. Removing it through the
  medium profile lost about 19% of interprotein CSMs and 34% of interprotein
  URPs, while providing no clear PPI-classification advantage. Retain the
  current `URPsupport` formulation in the large profile; do not claim that it
  solves PPI-level under-reporting.

## 2. Publish the stable MS2 workflow

- [x] Freeze the current validated MS2 defaults: automatic complexity profile,
  training-only Score.Diff prefilter when supported by the data, linear SVM,
  one hyperparameter search, three selected-model score repeats, and explicit
  post-scoring evidence polishing.
- [x] Complete and deliver the full Astral peptide-pair analysis. The report
  includes pre-classification and classified plots, xiView screenshots,
  MS-Viewer result keys, an Excel workbook, crosslinker comparisons, and the
  effects of instrument type and maximum-peak settings. DSSO performed more
  convincingly than DSBSO in the current Prospector workflow; the 500-peak
  searches likely suffered from Prospector's unmatched-peak penalty and do not
  resolve optimal Astral scoring parameters.
- [x] Rebuild the ribosome example as a concise, reproducible end-to-end
  demonstration of input, training, model inspection, result preparation,
  classification, polishing, and reporting.
- [x] Rewrite the GitHub README around the working ribosome example, with clear
  installation instructions, expected outputs, and links to more detailed
  documentation where appropriate.
- [ ] Run the README example from a clean R session and fresh package install,
  run the complete package checks, then tag the stable public MS2 release.

## 3. Improve PPI scoring and reporting after the MS2 release

- [x] Complete the bounded E. coli PPI adequacy audit using the 2025 Nature
  Methods resubmission, Kojak calls, human entrapment proteins, orthogonal
  SEC co-fractionation evidence, and STRING scores. The current Touchstone PPI
  list is well calibrated and strongly supported, but it under-recovers PPIs
  that Kojak reports and that independent evidence suggests are plausible.
- [x] Test whether stronger or absent `URPsupport` resolves PPI under-recovery.
  Stronger support reduced recovery, whereas omitting the feature greatly
  reduced CSM/URP recovery and did not clearly improve efficient PPI
  classification. Leave the CSM/URP feature unchanged.
- [ ] Define a dedicated one-row-per-PPI results table plus an associated
  evidence mapping that lists contributing URPs/CSMs without representing the
  PPI solely through its best CSM. Defer complex protein-inference resolution;
  retain transparent accession alternatives where feasible.
- [ ] Prototype a transparent secondary PPI evidence score based on the best
  URP score, additional independently positioned URPs, their score-weighted
  support, and CSM evidence. Continue training the primary SVM at CSM level;
  do not introduce separate SVMs for every summarization level without evidence
  that the simpler secondary aggregation fails.
- [ ] Validate any PPI score with decoys/entrapments, Kojak agreement,
  co-fractionation, and STRING. Treat co-fractionation and STRING strictly as
  held-out validation evidence rather than model features. Check portability on
  at least one additional dataset before changing the public PPI workflow.
- [ ] Revisit Bayesian integration and linear-peptide protein-ID priors only if
  the transparent aggregation score is inadequate.
- [ ] Deprioritized side quest: investigate high-confidence Kojak spectra for
  which Prospector reports no corresponding annotation, separately from cases
  where both engines identify the crosslink but Touchstone places it below the
  reporting threshold.

## 4. Repair and validate MS3 reconstruction

- [ ] Inventory the scan-linking and reconstruction code in `R/linkedScans.R`
  and document the assumptions made about scan relationships.
- [ ] Replace ad hoc scan matching with explicit validated joins and preserve
  reconstruction provenance in canonical flat columns.
- [ ] Add small fixtures covering missing scans, ambiguous links, charge/mass
  checks, and successful reconstruction.
- [ ] Keep reconstructed-MS3 support clearly marked experimental until those
  cases pass end-to-end tests.

## 5. Prospector integration and deferred extensions

- [ ] Improve integration between Protein Prospector Search Compare and
  Touchstone. Detect and canonicalize the optional MS-Product-derived fields
  for distinct product-ion cleavages, sequential/gapped ladders, and percent
  bond cleavage; report clearly when those annotations were not requested in
  Search Compare rather than treating their absence as failed evidence. Keep
  product-ion evidence restricted to the established backbone series and their
  cleavable-crosslinker variants (b/c and y/z with the supported `*`/`#`
  notation), collapsed to unique cleavage positions. Do not count neutral-loss
  variants, a ions, precursor/MH ions, or P-type diagnostic ions. Preserve the
  existing parser behavior until representative Prospector strings and
  regression fixtures establish any nomenclature change.
- [ ] Define reporting policies that degrade safely when MS-Product annotations
  are unavailable. Statistical classification must remain usable without these
  optional fields, while ion- or ladder-based polishing records whether it was
  applied, skipped for the entire dataset, or unavailable for particular rows.
  Allow an explicitly requested `Score.Diff >= 5` polishing fallback when ion
  annotations are unavailable; record the substitution rather than switching
  silently. Provide concise instructions for generating the required Search
  Compare output and revisit more direct Prospector-to-Touchstone transfer if a
  stable interface becomes available.
- [ ] Add experimental ternary-crosslink input only after the binary and MS3
  paths are stable.
- [ ] Define adapters from other search engines into Touchstone's canonical
  column names when comparison or rescoring work makes this worthwhile.
- [ ] Consider annotating each inter-protein pair with whether zero, one, or
  both constituent proteins have intra-protein crosslink support. This likely
  belongs in `calculatePairs()` as evidence annotation, not as an immediate
  reporting filter.
- [ ] Reconsider sequential ladder compliance as a possible model feature only
  if stringent manual polishing is repeatedly needed in independent datasets.
  Until then, keep it separate from statistical classification as a transparent
  evidence-quality policy rather than tuning the SVM to reproduce the TRiC
  reporting rule.
- [ ] Defer formal protein inference. In particular, do not attempt to resolve
  overlapping but non-identical accession sets until the stable CSM/URP
  workflow and PPI evidence-table requirements make the scope concrete.
- [ ] Decide whether to remove or clearly quarantine the deprecated
  `trainClassifier()` workflow and unused parallel (`furrr`/`future`) paths.
- [ ] Add continuous package checks and decide how to handle the large bundled
  example data that currently produces an `R CMD check` size note.
