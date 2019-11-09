## First submission

This is the first time the package is submitted to CRAN (second attempt).

## Test environments
* ubuntu xenial 16.04 (on travis-ci), r-release (3.6.1) and r-devel (4.0.0)
* mac osx xcode 9.4 (on travis-ci), r-release (3.6.1); could not build for r-devel due to https://cloud.r-project.org/bin/macosx/el-capitan/contrib/r-devel/ package repository being empty
* win-builder (in RStudio), r-release (3.6.1) and r-devel (4.0.0)

## R CMD check results
No ERRORs or WARNINGs from any checks.

* The Windows check returned one NOTE containing two possible issues:

1. 
> Possibly mis-spelled words in DESCRIPTION:  
> Heteroskedasticity (3:8)  
> Homoskedasticity (10:5)  
> heteroskedasticity (8:15, 11:49)  

This is not a misspelling. Alfredo R. Paloyo has [published an article](https://www.rwi-essen.de/media/content/pages/publikationen/ruhr-economic-papers/REP_11_300.pdf) on the etymology of the word, concluding that 'heteroskedasticity' 
is the correct spelling while 'heteroscedasticity' is incorrectly transliterated.

2. 
> Reading CITATION file fails with  
> there is no package called ‘skedastic’  
> when package is not installed.  

This is not a serious problem; in fact `citation("package")` appears to throw an 
error for most CRAN packages if the package is not installed.

* The Debian check on the previous submission attempt returned one additional NOTE:

> checking for detritus in the temp directory ... NOTE  
> Found the following files/directories:  
> ‘RtmpLz8nMX\hetplot’  
>  ‘RtmpLz8nMX\hetplot\X1_res_2_2019_11_07_11_52_44.png’  
> ...

These image files are written to a subdirectory of `tempdir()` by a function `hetplot` that creates heteroskedasticity plots. This is in line with the CRAN Repository Policy which states that   

> Packages should not write in the user’s home filespace (including clipboards), nor   
> anywhere else on the file system apart from the R session’s temporary directory

However, to avoid this issue recurring the relevant example has been changed to 
`## NOT RUN`
