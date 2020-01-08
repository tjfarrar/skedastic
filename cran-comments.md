## First submission

This is the first time the package is submitted to CRAN (fourth attempt).

## Revisions after previous submission
* Code was added to `hetplot` function to reset graphical parameters on exit

## Test environments
* ubuntu xenial 16.04 (on travis-ci), r-release (3.6.2) and r-devel (4.0.0)
* win-builder (in RStudio), r-release (3.6.2) and r-devel (4.0.0) 
* mac osx xcode 9.4 (on travis-ci), r-release (3.6.2); could not build for r-devel due to https://cloud.r-project.org/bin/macosx/el-capitan/contrib/r-devel/ package repository being empty

## R CMD check results
No ERRORs or WARNINGs from any checks.

* The Windows check returned one NOTE:

> Possibly mis-spelled words in DESCRIPTION:
>  BAMSET (14:69)
>  BLUS (16:9)
>  Breusch (18:40)
>  Glejser (13:8)
>  Goldfeld (12:5)
>  Heteroskedasticity (3:8)
>  Homoskedasticity (23:41)
>  Koenker (20:17)
>  Mittelhammer (14:5)
>  Quandt (12:18)
>  Theil (16:35)
>  Weisberg (21:73)
>  Yao (22:64)
>  heteroskedasticity (9:56, 25:21)

All of these are proper names or acronyms except for heteroskedasticity and 
homoskedasticity. The latter are not misspelled. Alfredo R. Paloyo has [published an article](https://www.rwi-essen.de/media/content/pages/publikationen/ruhr-economic-papers/REP_11_300.pdf) on the etymology of the word, concluding that 'heteroskedasticity' 
is the correct spelling while 'heteroscedasticity' is incorrectly transliterated.
