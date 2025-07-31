Touchstone CLMS re-scoring, classification, and more for Protein
Prospector
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# touchstone

<!-- badges: start -->
<!-- badges: end -->

The Touchstone library exists to improve the functionality of [Protein
Prospector](https://prospector.ucsf.edu/prospector/mshome.htm) for
Crosslinking Mass Spectrometry. Touchstone re-scores Prospector CLMS
results using an Support Vector Machine (SVM) classifier. It allows the
user to classify datasets at a desired False Discovery Rate (FDR)
threshold at various summarization levels: Crosslinked Spectral Matches
(CSMs), Unique Residue Pairs (URPs), or Protein Pairs (PPs).
Additionally Touchstone contains features to help with dataset
valdiation by measuring euclidean distances of crosslinks against
high-res structure files, or by querying
[STRING-db](https://string-db.org/) for String Scores of putative
protein interactions.

Touchstone can classify data into “Modules” which can either be used to
designate either domains within a larger polypeptide or stable
assemblies of multime polypeptides. There are convenience functions to
export data to [XiNet](https://crosslinkviewer.org/) and some internal
plotting functions for visualizing data at various summarization levels.

Touchstone is an R package that came out of my own needs when process a
lot of CLMS datasets searched with Prospector as a project scientist
supporting numerous projects over a 10-year period. It is therefore a
bit niche and wasn’t developed with a wide user base in mind.It is rough
around the edges and poorly supported. I am sharing it here because it
might be helpful to some users, but if you are looking for a smooth user
experience that doesn’t require tinkering in R, you might be better
served by other database search and re-scoring software.

I am currently only distributing a version here that is run in an R
command line environment, typically in RStudio. A demo version of a
graphical interface (built in Shiny) exists and is under development for
eventual integration with Protein Prospector. This will be more
accesible to a wider user-base.

## Installation

You can install the development version of touchstone from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mtrnka/Touchstone")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(touchstone)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
