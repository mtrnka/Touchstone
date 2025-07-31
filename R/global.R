the <- new.env(parent = emptyenv())

#default size of decoy database relative to target database
the$decoyScalingFactor = 1

#SVM model parameters:
params.best <- c("Score.Diff", "percMatched", "massError",
                 "z", "wtURP", "wtCSM", "xlinkClass",
                 "Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2")
params.noPercBond <- c("Score.Diff", "percMatched", "massError",
                 "z", "wtURP", "wtCSM", "xlinkClass")
params.best.nop <- c("Score.Diff", "percMatched", "massError",
                     "z", "wtCSM")

atomic.weight.da <- list(
# Values from:
# https://www.ciaaw.org/atomic-weights.htm
  "proton" = 1.0072765,
  "hydrogen" = 1.0078250,
  "oxygen" = 15.9949146,
  "sulfur" = 31.9720712,
  "nitrogen" = 14.0030740,
  "phosphorus" = 30.9737620
  )

H2O = 2 * atomic.weight.da[["proton"]] + atomic.weight.da[["oxygen"]]
SOH2 = H2O + atomic.weight.da[["sulfur"]]

# Path configuration (for shiny version):
exDir <- c(wd= '/mnt/pipeline/projects')
linkToXiView <- 'http://lanhuang2.physiology.uci.edu/crosslink-viewer/demo/Demo2.html'

# For generating MS-Viewer Links Correctly (mostly for shiny version):
queryTemplate.ms2 <- "
http://lanhuang2.physiology.uci.edu/prospector/cgi-bin/mssearch.cgi?
search_name=msproduct&
output_type=HTML&
report_title=MS-Product&version=6.2.29&
data_source=Data%20From%20File&
data_filename=%2Fmnt%2Fpipeline%2Fprojects%2FDSSOstar_rRibo%2Fsc%2FrRibo_DSSO_star_HCD%2FtstoneMS2.1%2FDSSOstar_rRibo%2FZ20200519-63_FTMSms2hcd.mgf&
use_instrument_ion_types=1&
msms_min_precursor_mass=0&
instrument_name=ESI-EThcD-high-res&display_graph=1&
msms_parent_mass_tolerance=10&
msms_parent_mass_tolerance_units=ppm&
fragment_masses_tolerance=30&
fragment_masses_tolerance_units=ppm&
msms_pk_filter=Max%20MSMS%20Pks&
msms_max_peaks=100&
scan_number=1000&
max_charge=4&
msms_precursor_charge=4&
sequence=SQK%28%2BDSG%29AIQDEIR&
s=1&
sequence2=Q%28Gln-%3Epyro-Glu%29QLPLPYEQLK%28%2BDSG%29HFYR&
s2=1&
count_pos_z=Ignore%20Basic%20AA&
link_search_type=No%20Link&
"

queryTemplate.ms2 <- unlist(stringr::str_split(stringr::str_replace_all(queryTemplate.ms2, "\\n", ""), "&"))
queryTemplate.ms2 <- stringr::str_split(queryTemplate.ms2, "=")
templateKeys.ms2 <- map_chr(queryTemplate.ms2, function(x) {x[1]})
templateVals.ms2 <- map_chr(queryTemplate.ms2, function(x) {x[2]})
templateVals.ms2 <- urltools::url_decode(templateVals.ms2)
names(templateVals.ms2) <- templateKeys.ms2
templateVals.ms2 <- templateVals.ms2[!is.na(templateVals.ms2)]

queryTemplate.ms3 <- "
http://lanhuang2.physiology.uci.edu/prospector/cgi-bin/mssearch.cgi?
search_name=msproduct&
output_type=HTML&
report_title=MS-Product&
version=6.2.29&
data_source=Data%20From%20File&
data_filename=
%2Fmnt%2Fpipeline%2Fprojects%2FDSSOms3_rRibo%2Fsc%2FrRibo_DSSO_ms3%2FtstoneMS3.1%2FDSSOms3_rRibo%2FZ20200519-49_ITMSms3cid.mgf&
use_instrument_ion_types=1&
msms_min_precursor_mass=0&
instrument_name=ESI-ION-TRAP-low-res&
display_graph=1&
msms_parent_mass_tolerance=10&
msms_parent_mass_tolerance_units=ppm&
fragment_masses_tolerance=0.7&
fragment_masses_tolerance_units=Da&
msms_pk_filter=Max%20MSMS%20Pks&
msms_max_peaks=40&
scan_number=37130&
max_charge=5&
msms_precursor_charge=5&
sequence=DHASIQMNVAEVDK%28XL:A-Alkene%29VTGR
count_pos_z=Ignore%20Basic%20AA&
s=1
"

queryTemplate.ms3 <- unlist(stringr::str_split(stringr::str_replace_all(queryTemplate.ms3, "\\n", ""), "&"))
queryTemplate.ms3 <- stringr::str_split(queryTemplate.ms3, "=")
templateKeys.ms3 <- map_chr(queryTemplate.ms3, function(x) {x[1]})
templateVals.ms3 <- map_chr(queryTemplate.ms3, function(x) {x[2]})
templateVals.ms3 <- urltools::url_decode(templateVals.ms3)
names(templateVals.ms3) <- templateKeys.ms3
templateVals.ms3 <- templateVals.ms3[!is.na(templateVals.ms3)]
