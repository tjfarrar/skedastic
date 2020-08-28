## Third submission

This is the third time the package is submitted to CRAN.

## Revisions after previous submission

* This is a very minor revision. There was a dependency on `plm` package, and I was notified that `plm` package has been earmarked for archiving. Thus I have removed the dependency.

## Test environments
* ubuntu xenial 16.04 (on travis-ci), r-release (4.0.2) and r-devel (4.1.0)
* win-builder (in RStudio), r-release (4.0.2) and r-devel (4.1.0) 
* mac osx xcode 9.4 (on travis-ci), r-release (4.0.2); could not build for r-devel due to https://cloud.r-project.org/bin/macosx/el-capitan/contrib/r-devel/ package repository being empty

## R CMD check results

* No ERRORs or WARNINGs from Windows (r-release, r-devel) checks or Mac OSX (r-release) check or Linux Xenial (r-release, r-devel) check.

* The Windows checks returned no NOTEs. Previous Windows checks have returned NOTEs about possible misspelled words or possible invalid DOIs but these NOTEs were spurious.
