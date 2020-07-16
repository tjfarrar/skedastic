## Second submission

This is the second time the package is submitted to CRAN.

## Revisions after previous submission
* More than 20 new functions were added to the package.

## Test environments
* ubuntu xenial 16.04 (on travis-ci), r-release (3.6.2) and r-devel (4.0.0)
* win-builder (in RStudio), r-release (3.6.2) and r-devel (4.0.0) 
* mac osx xcode 9.4 (on travis-ci), r-release (3.6.2); could not build for r-devel due to https://cloud.r-project.org/bin/macosx/el-capitan/contrib/r-devel/ package repository being empty

## R CMD check results
No ERRORs or WARNINGs from any checks.

* The Windows check returned a NOTE:

> Possibly mis-spelled words in DESCRIPTION:
>  Anscombe (11:42)
>  Bai (31:5)
>  Bickel (13:72)
>  Carapeto (16:68)
>  Cribari (34:20)
>  Cysneiros (34:5)
>  Diblasi (18:70)
>  Dufour (19:60)
>  Genest (20:18)
>  Godfrey (24:65)
>  Keselman (37:39)
>  Khalaf (19:68)
>  Kulinskaya (42:5)
>  Lehmann (43:38)
>  Neto (34:28)
>  Orme (25:5)
>  Rackauskas (31:63)
>  Simonoff (32:52)
>  Szroeter (35:5)
>  Tsai (32:65)
>  Verbyla (35:44)
>  Zhou (38:42)
>  Zuokas (32:5)

All of these are proper names of authors and not misspellings.

* The Windows check returned a NOTE:

> Found the following (possibly) invalid DOIs:
>  DOI: 10.2307/1911963
>    From: DESCRIPTION
>    Status: Forbidden
>    Message: 403
>  DOI: 10.2307/1912934
>    From: DESCRIPTION
>    Status: Forbidden
>    Message: 403
>  DOI: 10.2307/1913831
>    From: DESCRIPTION
>    Status: Forbidden
>    Message: 403
>  DOI: 10.2307/1913974
>    From: DESCRIPTION
>    Status: Forbidden
>    Message: 403
>  DOI: 10.2307/2986026
>    From: DESCRIPTION
>    Status: Forbidden
>    Message: 403

These DOIs have all been checked and are all valid: (10.2307/1911963)[https://doi.org/10.2307/1911963], (10.2307/1912934)[https://doi.org/10.2307/1912934], (10.2307/1913831)[https://doi.org/10.2307/1913831], 
(10.2307/1913974)[https://doi.org/10.2307/1913974],
(10.2307/2986026)[https://doi.org/10.2307/2986026].

