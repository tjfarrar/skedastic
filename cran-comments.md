## Fifth submission

This is the fourth time the package is submitted to CRAN.

## Revisions after previous submission

* The package is being resubmitted because it has a strong dependency on the `berryFunctions` package, which has been scheduled for archiving. This dependency has now been removed.

## Test environments

* ubuntu xenial 16.04 (on travis-ci), r-release and r-devel
* win-builder (in RStudio), r-oldrelease, r-release and r-devel

## R CMD check results

* No ERRORs or WARNINGs from Windows (r-release, r-devel, r-oldrelease) checks or Linux Xenial (r-release, r-devel) checks.

* The Windows checks returned one NOTE. This referred to possible misspelled words, but none of these words were misspelled. Most of them were surnames or acronyms, while 'heteroskedasticity' is a legitimate spelling of this word used in the academic literature (though sometimes spelled 'heteroscedasticity').

* The Windows checks also reported `Status: Forbidden`, `Message: 403` for some of the URLs and DOIs in the DESCRIPTION file. However, all of these URLs and DOIs are correct and working.
