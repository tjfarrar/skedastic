## First submission

This is the first time the package is submitted to CRAN (third attempt).

## Revisions after previous submission
* The copyright holder was added to the `Authors@R` field in DESCRIPTION
* References (with DOI or ISBN) were added to DESCRIPTION
* A data file `dpeakdat.Rda` was updated to a list of 500 objects instead of 300, 
  as a slow computation of 200 more objects had completed in the meantime.

## Test environments
* ubuntu xenial 16.04 (on travis-ci), r-release (3.6.1) and r-devel (4.0.0)
* mac osx xcode 9.4 (on travis-ci), r-release (3.6.1); could not build for r-devel due to https://cloud.r-project.org/bin/macosx/el-capitan/contrib/r-devel/ package repository being empty
* win-builder (in RStudio), r-release (3.6.1) and r-devel (4.0.0)

## R CMD check results
No ERRORs or WARNINGs from any checks.

* The Windows check returned one NOTE:

> Possibly mis-spelled words in DESCRIPTION:  
> Heteroskedasticity (3:8)  
> Homoskedasticity (10:5)  
> heteroskedasticity (8:15, 11:49)  

This is not a misspelling. Alfredo R. Paloyo has [published an article](https://www.rwi-essen.de/media/content/pages/publikationen/ruhr-economic-papers/REP_11_300.pdf) on the etymology of the word, concluding that 'heteroskedasticity' 
is the correct spelling while 'heteroscedasticity' is incorrectly transliterated.
