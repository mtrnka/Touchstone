# Touchstone stabilization roadmap

This roadmap records the agreed direction for making the current package a
stable successor to the original script-based Touchstone workflow. The public
data model will remain a flat data frame. More complicated internal structures
may be used where helpful, but package functions should return canonical flat
tables with enough provenance to understand how each result was produced.

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
  chosen profile. Automatic selection counts proteins with plausible repeated
  intra-protein CSM evidence rather than every reported accession, and retains
  the original dominant-protein safeguard for small systems with background
  matches.
- [x] Add explicit training and result objects, model-selection diagnostics,
  score-correlation plots, and reproducible decoy downsampling while retaining
  direct access to all fitted candidates.
- [ ] Explicitly evaluate a small set of prefilters, including `Score.Diff` and
  minimum product-ion evidence, using validation data rather than choosing the
  setting with the largest apparent yield. Score.Diff prefilter selection is
  now enabled for the large-system profile; product-ion filtering remains.
- [x] Make a linear SVM the default model. Radial candidates are available only
  for explicit, inspection-only experiments until independent validation shows
  a reproducible, material benefit.
- [ ] Complete the candidate audit. The candidate table now records each
  model's filters, features, parameters, validation FDR/yield, relative
  recovery within its kernel family, eligibility, recommendation status, and
  selection or rejection reason. Explicit repeat-split stability measurements
  remain to be designed and validated.
- [x] Separate statistical classification from evidence polishing. Provide
  named, transparent polishing policies (for example, minimum product ions or
  backbone-ladder coverage), record the applied rules and their before/after
  CSM counts, and recalculate FDR on the polished result.
- [x] Add the classified reporting path: apply the reporting-level thresholds
  to scored CSMs, rerun `calculatePairs()` on the surviving CSMs, and only then
  summarize to URP, peptide-pair, protein-pair, or module-pair. This ensures
  `numCSM` and `numURP` describe threshold-passing evidence; the training-only
  weighted count features are omitted from classified reporting tables.
- [ ] Validate the automated procedure on datasets spanning the three
  complexity profiles, including difficult DSSO data, before declaring the
  interface stable.

## 2. Repair and validate MS3 reconstruction

- [ ] Inventory the scan-linking and reconstruction code in `R/linkedScans.R`
  and document the assumptions made about scan relationships.
- [ ] Replace ad hoc scan matching with explicit validated joins and preserve
  reconstruction provenance in canonical flat columns.
- [ ] Add small fixtures covering missing scans, ambiguous links, charge/mass
  checks, and successful reconstruction.
- [ ] Keep reconstructed-MS3 support clearly marked experimental until those
  cases pass end-to-end tests.

## 3. Deferred extensions and cleanup

- [ ] Add experimental ternary-crosslink input only after the binary and MS3
  paths are stable.
- [ ] Define adapters from other search engines into Touchstone's canonical
  column names when comparison or rescoring work makes this worthwhile.
- [ ] Decide whether to remove or clearly quarantine the deprecated
  `trainClassifier()` workflow and unused parallel (`furrr`/`future`) paths.
- [ ] Refresh the README example after the stable training and polishing APIs
  exist.
- [ ] Add continuous package checks and decide how to handle the large bundled
  example data that currently produces an `R CMD check` size note.
