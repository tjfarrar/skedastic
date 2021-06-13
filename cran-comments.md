## Fourth submission

This is the fourth time the package is submitted to CRAN.

## Revisions after previous submission

* The package is being resubmitted because it was archived on 2021-06-08 due to an issue not having been corrected quickly enough. The issue was that the package had a dependency on `het.test` package which was archived. This dependency has now been removed.

## Test environments

* ubuntu xenial 16.04 (on travis-ci), r-release and r-devel
* win-builder (in RStudio), r-oldrelease, r-release and r-devel 
* mac osx xcode 9.4 (on travis-ci), r-release ; could not build for r-devel due to https://cloud.r-project.org/bin/macosx/el-capitan/contrib/r-devel/ package repository being empty

## R CMD check results

* No ERRORs or WARNINGs from Windows (r-release, r-devel) checks or Mac OSX (r-release) check or Linux Xenial (r-release, r-devel) check.

* The Windows checks returned one NOTE. This referred to possible misspelled words, but none of these words were misspelled. Most of them were surnames or acronyms, while 'heteroskedasticity' is a legitimate spelling of this word used in the academic literature (though sometimes spelled 'heteroscedasticity').
