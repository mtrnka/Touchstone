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
- [ ] Explicitly evaluate a small set of prefilters, including `Score.Diff` and
  minimum product-ion evidence, using validation data rather than choosing the
  setting with the largest apparent yield. Score.Diff prefilter selection is
  is evaluated whenever there are enough target and decoy CSMs for comparison,
  independently of the feature-complexity profile. It limits model fitting
  without removing lower-Score.Diff CSMs from subsequent scoring. Product-ion
  evidence remains a reporting-polish option rather than a training prefilter.
  The legacy low-FDR interprotein summary (`interInt`) is retained internally
  for this prefilter choice pending validation on the original E. coli data,
  but is no longer exposed in the hyperparameter candidate audit.
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
- [ ] Re-test the training-only Score.Diff prefilter on the large E. coli data
  used to develop the original procedure. Confirm that the automatic choice can
  recover the historically useful 15--20 training thresholds and that scoring
  the complete input recovers credible lower-Score.Diff CSMs without degrading
  FDR or manual quality.
- [x] Replace the self-inclusive, hard-thresholded `wtCSM` and `wtURP` training
  features with `CSMsupport` and `URPsupport`. The default complexity profiles
  now use self-excluded corroborating evidence with a bounded soft contribution.
  The legacy features remain available for compatibility and comparison but are
  omitted from the new automatic profiles and classified reporting tables.

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

## 4. Deferred extensions and cleanup

- [ ] Add experimental ternary-crosslink input only after the binary and MS3
  paths are stable.
- [ ] Define adapters from other search engines into Touchstone's canonical
  column names when comparison or rescoring work makes this worthwhile.
- [ ] Consider annotating each inter-protein pair with whether zero, one, or
  both constituent proteins have intra-protein crosslink support. This likely
  belongs in `calculatePairs()` as evidence annotation, not as an immediate
  reporting filter.
- [ ] Decide whether to remove or clearly quarantine the deprecated
  `trainClassifier()` workflow and unused parallel (`furrr`/`future`) paths.
- [ ] Add continuous package checks and decide how to handle the large bundled
  example data that currently produces an `R CMD check` size note.
