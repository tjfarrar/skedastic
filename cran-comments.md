## Sixth submission

This is the sixth time the package is submitted to CRAN.

## Revisions after previous submission

* A substantial amount of new functionality has been added to the package, particularly through the addition of the functions `alvm.fit` and `anlvm.fit`.

## Test environments

* ubuntu xenial 16.04 (on travis-ci), r-release and r-devel
* win-builder (in RStudio), r-oldrelease, r-release and r-devel

## R CMD check results

* No ERRORs or WARNINGs from Windows (r-release, r-devel, r-oldrelease) checks
* Linux Xenial (r-release, r-devel) checks produced ERRORs, but these were all related to the failture of packages such as `devtools` and `pkgdown` to load, and therefore not the result of problems with the `skedastic` package

* All of the Windows checks returned one NOTE. This referred to `Status: Forbidden`, `Message: 403` for some of the URLs and DOIs in the DESCRIPTION file. However, all of these URLs and DOIs have been checked manually and are correct. It also referred to possible spelling errors in the DESCRIPTION file, but these were all either proper names, Latin expressions like et al., or correctly spelt technical terms like heteroskedasticity.

* The r-devel check returned a second NOTE, `checkRd: (-1) GSS.Rd:86: Escaped LaTeX specials: \&`. However, the `GSS.Rd` file does not contain this character.
