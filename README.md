Touchstone: CLMS re-scoring, classification, biological inference and
more for Protein Prospector crosslink searches
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Protein Prospector - Touchstone CLMS pipeline

<!-- badges: start -->
<!-- badges: end -->

The Touchstone library exists to improve the functionality of
<a href="https://prospector.ucsf.edu/prospector/mshome.htm"
target="_blank">Protein Prospector</a> for Crosslinking Mass
Spectrometry (Trnka et al. 2014). Touchstone re-scores Prospector CLMS
results using an Support Vector Machine (SVM) classifier that does a
better job discriminating between correct and incorrect crosslinks than
the Prospector scores.

It allows the user to classify datasets at a desired False Discovery
Rate (FDR) threshold at various summarization levels: Crosslinked
Spectral Matches (CSMs), Unique Residue Pairs (URPs), or Protein Pairs
(PPs). Touchstone’s FDR assessment are highly consistent with error
rates measured assessed by various benchmarking datasets (Beveridge et
al. 2020; Matzinger et al. 2022; Fischer et al. 2025).

Additionally Touchstone contains features to help with dataset
valdiation by measuring euclidean distances of crosslinks against
high-res structure files, or by querying
<a href="https://string-db.org/" target="_blank">STRING-db</a> for
String Scores of putative protein interactions.

To aid with biological inference, Touchstone can optionally classify
data into “Modules” which can designate either domains within a larger
polypeptide or stable assemblies of multiple polypeptides (or both).
There are convenience functions to export data to
<a href="https://crosslinkviewer.org/" target="_blank">XiNet</a> and
some internal plotting functions for quantitating CSMs across proteins
or modules.

Touchstone is an R package that that I developed to address my own needs
when analyzing CLMS datasets searched with Prospector, as a project
scientist supporting numerous <a
href="https://scholar.google.com/citations?hl=en&amp;user=Gae1r_AAAAAJ&amp;view_op=list_works&amp;sortby=pubdate"
target="_blank">projects</a> over the last decade or so. It is therefore
a bit niche and wasn’t developed with a wide user base in mind. Nor am I
a software developer. So, it is rough around the edges. I am sharing it
here because it might be helpful to some users, but if you are looking
for a smooth user experience that doesn’t require tinkering in R, you
might be better served by other CLMS database search and re-scoring
software.

I am currently only distributing a version here that is run in an R
command line environment, typically in
<a href="https://posit.co/download/rstudio-desktop/"
target="_blank">RStudio</a>. A demo version of a graphical interface
(built in <a href="https://shiny.posit.co/" target="_blank">Shiny</a>)
exists and is under development for eventual integration with Protein
Prospector. This will be more accesible to a wider user-base. The
eventual goal is to make Prospector CLMS searches more widely accessible
to the research community.

The graphical demo version can be accessed at
<a href="https://prospts.shinyapps.io/tstoneapp/"
target="_blank">shinyapps.io</a>

Additional instructions for the graphical demo are <a
href="https://msf.ucsf.edu/mike/crosslinkingClass/dataset_summary.html"
target="_blank">here</a>

The rest of this mini-vignette will refer to running Touchstone inside
of RStudio.

## Installation

The touchstone library is distributed on
<a href="https://github.com/" target="_blank">GitHub</a>. Install with:

``` r
# install.packages("devtools")
devtools::install_github("mtrnka/Touchstone")
```

## 80S Ribosome data acquired by MS2.HCD.

80S ribosome was produced using a rabbit reticulocyte cell free
expression system. 80S ribosomes were crosslinked with the cleavable
reagent DSSO (Kao et al. 2011).

This ribosomal system for method development and optimization of CLMS
workflows. There is a high-res EM structure of the complex which can be
helpful in determining if the crosslinked are assigned correctly or not,
<a href="https://www.rcsb.org/structure/6hcj"
target="_blank">pdb:6HCJ</a>.

The example dataset included with Touchstone is from 80S sample analyzed
using a stepped-HCD MS2 acquisition cycle (dataset \#1). The modulefile
categorizes the 80 or so ribosomal proteins to either the large (60S) or
small (40S) subunits.

<figure>
<img src="https://cdn.rcsb.org/images/structures/6hcj_assembly-1.jpeg"
style="width:30.0%" alt="cryoEM structure of Rabbit 80S ribosome" />
<figcaption aria-hidden="true">cryoEM structure of Rabbit 80S
ribosome</figcaption>
</figure>

``` r
library(touchstone)
library(tidyverse)

pathToDemoFile <- system.file("extdata", "rRibo_DSSO_sthcd_scOut.txt", package="touchstone")
ribo <- readProspectorXLOutput(pathToDemoFile, minPepLen = 4, minIons = 0)
```

The crosslinking database search was run with Prospector program **Batch
Tag**. **Search Compare** is then used to process the search results and
create a tab delimited text file compatible with Touchstone processing.
Touchstone expects unclassified search results with certain parameters
included in the report. The recommended Search Compare Paramters are
inlcuded as an example file:

``` r

# install.packages("jsonlite")
touchstone_example("tstoneMS2.4.json") %>% jsonlite::read_json()
#> $creation_time
#> [1] "Fri Aug 20 13:50:37 2021"
#> 
#> $pp_version
#> [1] "6.3.23"
#> 
#> $parameters
#> $parameters$output_directory
#> [1] "/home/mtrnka/rRibo/sc"
#> 
#> $parameters$output_filename
#> [1] "__outputfilename__"
#> 
#> $parameters$results_full_paths
#> [1] "__fullpathtobatchtagresults__"
#> 
#> $parameters$results_names
#> [1] "doesItMatter"
#> 
#> $parameters$project_names
#> [1] "__fplfolder__"
#> 
#> $parameters$area_threshold
#> [1] "0"
#> 
#> $parameters$best_disc_only
#> [1] "1"
#> 
#> $parameters$cache_name
#> [1] "results"
#> 
#> $parameters$comp_mask_type
#> [1] "OR"
#> 
#> $parameters$disc_score_graph
#> [1] "1"
#> 
#> $parameters$intensity_threshold
#> [1] "0"
#> 
#> $parameters$library_options
#> [1] "Remove%20Ambiguous%20Peptides"
#> 
#> $parameters$matched_intensity
#> [1] "1"
#> 
#> $parameters$base_intensity
#> [1] "1"
#> 
#> $parameters$msms_err
#> [1] "1"
#> 
#> $parameters$msms_int
#> [1] "1"
#> 
#> $parameters$msms_ions
#> [1] "1"
#> 
#> $parameters$msms_mz
#> [1] "1"
#> 
#> $parameters$max_peptide_evalue
#> [1] "1000"
#> 
#> $parameters$max_peptide_fdr
#> [1] "1"
#> 
#> $parameters$max_protein_evalue
#> [1] "1000"
#> 
#> $parameters$max_protein_fdr
#> [1] "5"
#> 
#> $parameters$merge_option
#> [1] "Separated"
#> 
#> $parameters$min_best_disc_score_ESI_ETD_high_res
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_ESI_ETD_low_res
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_ESI_EThcD_high_res
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_ESI_FT_ICR_CID
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_ESI_FT_ICR_ECD
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_ESI_ION_TRAP_low_res
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_ESI_Q_TOF
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_ESI_Q_high_res
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_MALDI_Q_TOF
#> [1] "0.0"
#> 
#> $parameters$min_best_disc_score_MALDI_TOFTOF
#> [1] "0.0"
#> 
#> $parameters$min_peptide_score
#> [1] "0"
#> 
#> $parameters$min_protein_score
#> [1] "5"
#> 
#> $parameters$msms_max_peaks
#> [1] "85"
#> 
#> $parameters$msms_pk_filter
#> [1] "Max%20MSMS%20Pks"
#> 
#> $parameters$orbitrap
#> [1] "1"
#> 
#> $parameters$output_filename
#> [1] "name"
#> 
#> $parameters$parent_mass_convert
#> [1] "monoisotopic"
#> 
#> $parameters$peak_list_type
#> [1] "mzML"
#> 
#> $parameters$peptide_fdr_type
#> [1] "Peptide"
#> 
#> $parameters$peptide_filter
#> [1] "Keep%20Replicates"
#> 
#> $parameters$peptide_mod_type
#> [1] "Variable%20Mods%20Only"
#> 
#> $parameters$pepxml_options
#> [1] "Spectral%20Top%20Hit"
#> 
#> $parameters$percent_C13
#> [1] "100"
#> 
#> $parameters$percent_H2
#> [1] "100"
#> 
#> $parameters$percent_N15
#> [1] "100"
#> 
#> $parameters$percent_O18
#> [1] "100"
#> 
#> $parameters$percent_bond_cleavage
#> [1] "1"
#> 
#> $parameters$purity_correction
#> [1] "From%20formulae"
#> 
#> $parameters$quan_batch_option
#> [1] "Write"
#> 
#> $parameters$quan_type
#> [1] "DTT_C%202H%20%28C%29"
#> 
#> $parameters$raw_type
#> [1] "MS%20Precursor"
#> 
#> $parameters$rep_q_n_sdv
#> [1] "2.0"
#> 
#> $parameters$report_accession
#> [1] "1"
#> 
#> $parameters$report_charge
#> [1] "1"
#> 
#> $parameters$report_db_peptide
#> [1] "1"
#> 
#> $parameters$report_elem_comp
#> [1] "1"
#> 
#> $parameters$report_end_aa
#> [1] "1"
#> 
#> $parameters$report_error
#> [1] "1"
#> 
#> $parameters$report_expectation
#> [1] "1"
#> 
#> $parameters$report_hits_type
#> [1] "Union"
#> 
#> $parameters$report_homologous_proteins
#> [1] "All"
#> 
#> $parameters$report_links
#> [1] "1"
#> 
#> $parameters$report_m_over_z
#> [1] "1"
#> 
#> $parameters$report_msms_info
#> [1] "1"
#> 
#> $parameters$report_name
#> [1] "1"
#> 
#> $parameters$report_next_aa
#> [1] "0"
#> 
#> $parameters$report_num_pks
#> [1] "1"
#> 
#> $parameters$report_number
#> [1] "1"
#> 
#> $parameters$report_previous_aa
#> [1] "0"
#> 
#> $parameters$report_prot_len
#> [1] "1"
#> 
#> $parameters$report_protein_mod
#> [1] "1"
#> 
#> $parameters$report_repeats
#> [1] "1"
#> 
#> $parameters$report_score
#> [1] "1"
#> 
#> $parameters$report_score_diff
#> [1] "1"
#> 
#> $parameters$report_species
#> [1] "1"
#> 
#> $parameters$report_start_aa
#> [1] "1"
#> 
#> $parameters$report_time
#> [1] "1"
#> 
#> $parameters$report_type
#> [1] "Crosslinked%20Peptides"
#> 
#> $parameters$report_unmatched
#> [1] "1"
#> 
#> $parameters$report_xl_aa
#> [1] "1"
#> 
#> $parameters$report_xl_expectation
#> [1] "1"
#> 
#> $parameters$report_xl_peptide
#> [1] "1"
#> 
#> $parameters$report_xl_rank
#> [1] "1"
#> 
#> $parameters$report_xl_score
#> [1] "1"
#> 
#> $parameters$reporter_ion_window
#> [1] "0.4"
#> 
#> $parameters$resolution
#> [1] "60000"
#> 
#> $parameters$results_file
#> [1] ""
#> 
#> $parameters$rt_int_end
#> [1] "30.0"
#> 
#> $parameters$rt_int_start
#> [1] "-10.0"
#> 
#> $parameters$save_format
#> [1] "MS-Viewer%20Files"
#> 
#> $parameters$score_limit_type
#> [1] "Score%20Limits%20Only"
#> 
#> $parameters$search_name
#> [1] "searchcompare"
#> 
#> $parameters$slip_threshold
#> [1] "6"
#> 
#> $parameters$snr_threshold
#> [1] "10.0"
#> 
#> $parameters$sort_type
#> [1] "Crosslink%20AA"
#> 
#> $parameters$uncompressed_peak_lists
#> [1] "1"
#> 
#> $parameters$version
#> [1] "6.3.23"
#> 
#> $parameters$write_cache
#> [1] "No%20Cache%20File"
#> 
#> $parameters$xl_max_low_expectation
#> [1] "1000"
#> 
#> $parameters$xl_min_low_score
#> [1] "0.0"
#> 
#> $parameters$xl_min_score_diff
#> [1] "0.0"
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

## References:

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-beveridge_synthetic_2020" class="csl-entry">

Beveridge, Rebecca, Johannes Stadlmann, Josef M. Penninger, and Karl
Mechtler. 2020. “A Synthetic Peptide Library for Benchmarking
Crosslinking-Mass Spectrometry Search Engines for Proteins and Protein
Complexes.” *Nature Communications* 11 (1): 742.
<https://doi.org/10.1038/s41467-020-14608-2>.

</div>

<div id="ref-fischer_assessment_2025" class="csl-entry">

Fischer, Lutz, Lars Kolbowski, Swantje Lenz, James E. Bruce, Robert J.
Chalkley, Michael R. Hoopmann, David D. Shteynberg, et al. 2025.
“Assessment of Reported Error Rates in Crosslinking Mass Spectrometry.”
bioRxiv. <https://doi.org/10.1101/2025.04.27.649519>.

</div>

<div id="ref-kao_development_2011" class="csl-entry">

Kao, Athit, Chi-li Chiu, Danielle Vellucci, Yingying Yang, Vishal R.
Patel, Shenheng Guan, Arlo Randall, Pierre Baldi, Scott D. Rychnovsky,
and Lan Huang. 2011. “Development of a Novel Cross-Linking Strategy for
Fast and Accurate Identification of Cross-Linked Peptides of Protein
Complexes \*.” *Molecular & Cellular Proteomics* 10 (1).
<https://doi.org/10.1074/mcp.M110.002212>.

</div>

<div id="ref-matzinger_mimicked_2022" class="csl-entry">

Matzinger, Manuel, Adrian Vasiu, Mathias Madalinski, Fränze Müller,
Florian Stanek, and Karl Mechtler. 2022. “Mimicked Synthetic Ribosomal
Protein Complex for Benchmarking Crosslinking Mass Spectrometry
Workflows.” *Nature Communications* 13 (1): 3975.
<https://doi.org/10.1038/s41467-022-31701-w>.

</div>

<div id="ref-trnka_matching_2014" class="csl-entry">

Trnka, Michael J., Peter R. Baker, Philip J. J. Robinson, A. L.
Burlingame, and Robert J. Chalkley. 2014. “Matching Cross-Linked Peptide
Spectra: Only as Good as the Worse Identification \*.” *Molecular &
Cellular Proteomics* 13 (2): 420–34.
<https://doi.org/10.1074/mcp.M113.034009>.

</div>

</div>
