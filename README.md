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
devtools::install_github("mtrnka/Touchstone@main")
#> Skipping install of 'touchstone' from a github remote, the SHA1 (9edc6295) has not changed since last install.
#>   Use `force = TRUE` to force installation
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

pathToDemoFile <- touchstone_example("rRibo_DSSO_sthcd_scOut.txt")
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
touchstone_example("tstoneMS2.4.json") %>% 
  jsonlite::read_json()
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
