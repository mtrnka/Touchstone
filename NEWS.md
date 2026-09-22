# touchstone 0.1.0

Initial stable release of the binary-crosslink MS2 workflow.

- Added reproducible `trainCrosslinkScore()` model selection with conservative
  linear-SVM defaults, grouped cross-validation, training-only Score.Diff
  prefiltering, dataset-complexity feature profiles, and a three-fit selected
  model score ensemble.
- Added explicit training and result objects with candidate audits,
  score-correlation diagnostics, stable threshold estimation, and retained
  access to fitted alternatives.
- Separated reporting-level preparation, statistical classification, and
  optional evidence polishing for CSM, unique-residue-pair, peptide-pair,
  protein-pair, and module-pair results.
- Added transparent product-ion and ladder polishing audits, classifier-aware
  exports, improved FDR plots, structural reporting, xiView/xiNet output,
  ChimeraX pseudobonds, and optional STRING evidence retrieval.
- Expanded regression tests for FDR calculations, pair construction,
  summarization, input-column normalization, model training, cross-fitting,
  plotting, result preparation, and exports.
- Rebuilt the README around a tested end-to-end rabbit ribosome analysis.
