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
>  BAMSET (13:69)
>  BLUS (15:9)
>  Breusch (17:40)
>  Glejser (12:8)
>  Goldfeld (11:5)
>  Heteroskedasticity (3:8)
>  Homoskedasticity (22:41)
>  Koenker (19:17)
>  Mittelhammer (13:5)
>  Quandt (11:18)
>  Theil (15:35)
>  Weisberg (20:73)
>  Yao (21:64)
>  heteroskedasticity (8:56, 24:21)

All of these are proper names or acronyms except for heteroskedasticity and 
homoskedasticity. The latter are not misspelled. Alfredo R. Paloyo has [published an article](https://www.rwi-essen.de/media/content/pages/publikationen/ruhr-economic-papers/REP_11_300.pdf) on the etymology of the word, concluding that 'heteroskedasticity' 
is the correct spelling while 'heteroscedasticity' is incorrectly transliterated.
