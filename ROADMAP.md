# Touchstone stabilization roadmap

This roadmap records the agreed direction for making the current package a
stable successor to the original script-based Touchstone workflow. The public
data model will remain a flat data frame. More complicated internal structures
may be used where helpful, but package functions should return canonical flat
tables with enough provenance to understand how each result was produced.

The current release priority is to complete section 1, evaluate the finalized
MS2 scoring workflow on the full Astral dataset, and report that analysis.
The ribosome example and GitHub README will then be completed before work begins
on MS3 reconstruction or the other deferred extensions.

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
  remains linear until independent validation shows a reproducible benefit.
- [ ] Complete the candidate audit. The candidate table now records each
  model's filters, features, parameters, validation FDR/yield, relative
  recovery within its kernel family, concise intra- and interprotein rank
  correlations over both the full range and high-Score.Diff tail, eligibility,
  recommendation status, and selection or rejection reason. Models with a
  weakest correlation below the kernel-specific minimum (0.2 for linear and
  0.5 for radial by default) are excluded before the near-best recovery range
  is calculated. Selection falls back to intraprotein recovery when credible
  models recover no interprotein links.
  Explicit repeat-split stability measurements remain to be designed and
  validated across the Astral, E. coli, and small-system test datasets.
- [x] Separate statistical classification from evidence polishing. Provide
  named, transparent polishing policies (for example, minimum product ions or
  backbone-ladder coverage), record the applied rules and their before/after
  CSM counts, and recalculate FDR on the polished result.
- [x] Add the classified reporting path: apply the reporting-level thresholds
  to scored CSMs, rerun `calculatePairs()` on the surviving CSMs, and only then
  summarize to URP, peptide-pair, protein-pair, or module-pair. This ensures
  `numCSM` and `numURP` describe threshold-passing evidence; the training-only
  weighted count features are omitted from classified reporting tables.
- [ ] Define and run an MS2 testing protocol on representative datasets spanning
  the small, medium, and large complexity profiles, including difficult DSSO
  data. Record the expected outputs, selected settings, model diagnostics,
  thresholds, recovery, and seed-to-seed stability before declaring the
  interface stable. The initial local harness is
  `scratch/touchstoneComplexityTest.R`; this item remains open until its dataset
  suite and repeated-seed comparisons have been reviewed.
- [x] Validate the current small-system workflow against the collaborator-
  delivered TRiC sample-1 result. The exact historical thresholding and ladder-
  polishing order reproduced the same total of 554 target URPs, sharing 277 of
  278 inter-protein and 272 of 276 intra-protein URPs. Retain the historical
  sequential-ladder requirement as an explicit conservative reporting policy;
  do not treat the high open/closed-structure violation rates as a useful truth
  set for further tuning.
- [ ] Complete validation of the HCD and EThcD translocon results. The initial
  4UG0 comparison found a higher mapped inter-protein violation rate for EThcD
  than HCD, so the larger EThcD result cannot yet be accepted solely as improved
  fragmentation. Trace the historical Keenan/eLife analysis, compare shared and
  fragmentation-method-specific links, and test whether EThcD-only recovery is
  supported by improved peptide-2 fragmentation or reflects a liberal score
  threshold. The same comparison supports de-emphasizing radial models, whose
  unique additions were more violation-prone; radial remains an explicitly
  requested exploratory option rather than a public default.
- [x] Re-test the training-only Score.Diff prefilter on the large E. coli data
  used to develop the original procedure. Fixed linear models at Score.Diff
  thresholds 0, 5, 10, 15, and 20 were scored against the complete input, with
  the historical minimum-four-product-ions-per-peptide rule applied before
  independent 2% CSM and protein-pair thresholds. Thresholds 5--20 produced a
  stable plateau. The automatic value of 10 recovered 10,437 normal-target
  inter-protein CSMs, 136 normal decoys, 310 total entrapment hits, and 342
  inter-protein PPIs, closely matching the historical 10,706/137/266 and about
  366 PPI result. Lower-Score.Diff CSM recovery was retained and reported.
- [x] Replace the self-inclusive, hard-thresholded `wtCSM` and `wtURP` training
  features with `CSMsupport` and `URPsupport`. The default complexity profiles
  now use self-excluded corroborating evidence with a bounded soft contribution.
  The legacy features remain available for compatibility and comparison but are
  omitted from the new automatic profiles and classified reporting tables.
- [ ] Perform a bounded E. coli PPI adequacy audit before changing the scoring
  architecture. Using the current accession-pair definitions, tabulate every
  candidate inter-protein PPI with its SVM score/FDR band, CSM and URP support,
  decoy/entrapment status, historical Touchstone and Kojak classification,
  co-fractionation evidence, and STRING score. Evaluate two separate questions:
  whether 2% PPI FDR is calibrated by entrapment, and whether orthogonal support
  declines sensibly below the Touchstone threshold. Treat co-fractionation and
  STRING as validation evidence, not training features.
- [ ] Use the PPI audit as a decision gate. Retain one CSM-level SVM score with
  independent CSM, URP, and PPI thresholds if its PPI ranking is adequate. If
  ranking is inadequate, first test a transparent linear PPI aggregation score
  based on the strongest and additional independent URPs. Do not introduce
  separate SVMs at every summarization level or a Bayesian model without
  evidence that the simpler design fails.

## 2. Analyze and present the stable MS2 workflow

- [ ] Freeze the validated MS2 defaults and use that version for the full
  Astral analysis.
- [ ] Complete and report the full Astral analysis before revising the public
  example.
- [ ] Rebuild the ribosome example as a concise, reproducible end-to-end
  demonstration of input, training, model inspection, result preparation,
  classification, polishing, and reporting.
- [ ] Rewrite the GitHub README around the working ribosome example, with clear
  installation instructions, expected outputs, and links to more detailed
  documentation where appropriate.
- [ ] Run the README example from a clean R session and fresh package install,
  then tag the stable public MS2 release.

## 3. Repair and validate MS3 reconstruction

- [ ] Inventory the scan-linking and reconstruction code in `R/linkedScans.R`
  and document the assumptions made about scan relationships.
- [ ] Replace ad hoc scan matching with explicit validated joins and preserve
  reconstruction provenance in canonical flat columns.
- [ ] Add small fixtures covering missing scans, ambiguous links, charge/mass
  checks, and successful reconstruction.
- [ ] Keep reconstructed-MS3 support clearly marked experimental until those
  cases pass end-to-end tests.

## 4. Prospector integration and deferred extensions

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
- [ ] Revisit a dedicated PPI results format only after the bounded PPI audit.
  A future format may separate one-row-per-PPI summaries from URP/CSM evidence
  mappings, but the current best-row representation remains sufficient for the
  stabilization analysis.
- [ ] Defer Bayesian evidence integration, linear-peptide priors, and formal
  protein inference. In particular, do not attempt to resolve overlapping but
  non-identical accession sets until the stable CSM/URP workflow and the need
  for a separate PPI score have been established.
- [ ] Decide whether to remove or clearly quarantine the deprecated
  `trainClassifier()` workflow and unused parallel (`furrr`/`future`) paths.
- [ ] Add continuous package checks and decide how to handle the large bundled
  example data that currently produces an `R CMD check` size note.
