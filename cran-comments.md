## Minor revision

The only change in this version of the package is that the `CITATION` file has been updated, as I forgot to update it when a substantially modified package version was published in 2022.

## Test environments

* ubuntu xenial 16.04 (on travis-ci), r-release and r-devel
* win-builder (in RStudio), r-oldrelease, r-release and r-devel

## R CMD check results

* No ERRORs or WARNINGs from Windows (r-release, r-devel, r-oldrelease) checks
* Linux Xenial (r-release, r-devel) checks produced ERRORs, but these were all related to the failture of packages such as `devtools` and `pkgdown` to load, and therefore not the result of problems with the `skedastic` package

* All of the Windows checks returned one NOTE. This referred to `Status: 500`, `Message: Internal Server Error` for one of the URLs in the DESCRIPTION file. However, the URL was checked manually in a web browser and it works.
