## Third submission

This is the third time the package is submitted to CRAN.

## Revisions after previous submission

* This is a very minor revision. There was a dependency on `plm` package, and I was notified that `plm` package has been earmarked for archiving. Thus I have removed the dependency.

## Test environments
* ubuntu xenial 16.04 (on travis-ci), r-release (4.0.2) and r-devel (4.1.0)
* win-builder (in RStudio), r-release (4.0.2) and r-devel (4.1.0) 
* mac osx xcode 9.4 (on travis-ci), r-release (4.0.2); could not build for r-devel due to https://cloud.r-project.org/bin/macosx/el-capitan/contrib/r-devel/ package repository being empty

## R CMD check results

* No ERRORs or WARNINGs from Windows (r-release, r-devel) checks or Mac OSX (r-release) check or Linux Xenial (r-devel) check. The Linux Xenial (r-release) check on travis-ci has an ERROR due to failure of lazy loading for some packages. This seems to be a temporary and unavoidable issue external to my package that will hopefully resolve itself.

* The Windows check returned a NOTE:

> Possibly mis-spelled words in DESCRIPTION:
>  Anscombe (11:42)
>  BAMSET (13:5)
>  BLUS (40:43)
>  Bai (31:5)
>  Bickel (13:72)
>  Breusch (14:42)
>  Carapeto (16:68)
>  Cribari (34:20)
>  Cysneiros (34:5)
>  Diblasi (18:70)
>  Dufour (19:60)
>  Genest (20:18)
>  Glejser (22:41)
>  Godfrey (24:65)
>  Goldfeld (25:50)
>  Heteroskedasticity (3:8)
>  Homoskedasticity (44:5)
>  Keselman (37:39)
>  Khalaf (19:68)
>  Koenker (16:17)
>  Kulinskaya (42:5)
>  Lehmann (43:38)
>  Mittelhammer (24:5)
>  Neto (34:28)
>  Orme (25:5)
>  Quandt (25:63)
>  Rackauskas (31:63)
>  Simonoff (32:52)
>  Szroeter (35:5)
>  Theil (40:61)
>  Tsai (32:65)
>  Verbyla (35:44)
>  Weisberg (17:61)
>  Yao (29:45)
>  Zhou (38:42)
>  Zuokas (32:5)
>  heteroskedasticity (9:56, 39:44, 45:49)

Most of these are proper names of authors and not misspellings. 'Heteroskedasticity' and 'homoskedasticity' are also not mispelled.

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

