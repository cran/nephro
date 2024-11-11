# Nephro with CKiD U25 equations
Xavier Gaeta MD, PhD

## Nephro - CKiD U25

This is fork of the R **nephro** package that’s hosted on CRAN.org which
includes functions to calculate pediatric eGFR based on the following
publication:

Pierce CB, Muñoz A, Ng DK, Warady BA, Furth SL, Schwartz GJ. **Age- and
sex-dependent clinical equations to estimate glomerular filtration rates
in children and young adults with chronic kidney disease.** *Kidney
International*. 2021;99(4):948–956.
[doi:10.1016/j.kint.2020.10.047](https://doi.org/10.1016/j.kint.2020.10.047 "Persistent DOI link for Pierce et al. paper")

## Installation

You can install the package from Github

``` r
devtools::install_github("xgaeta/nephro-CKiD")
```

then include the package in your script:

``` r
library(nephro)
```

## Using the functions

These new equations follow the same style as the rest of the package,
including sex as a binary variable with 0 for female and 1 for male.
Height is provided in cm. They currently use US based equations. Outputs
are in mL/min/1.73m<sup>2</sup>.

``` r
CKiD.U25.cystatin(cystatin = 1.2, age = 9.5, sex = 0)
```

    [1] 68.4

``` r
CKiD.U25.creatinine(creatinine = 0.8, age = 9.5, sex = 0, ht = 132)
```

    [1] 63.1

``` r
CKiD.U25.combined(creatinine = 0.8, cystatin = 1.2, age = 9.5, sex = 0, ht = 132, verbose = TRUE)
```

         eGFRU25.cr eGFRU25.cys sGFRU25.avg
    [1,]       63.1        68.4        65.8

Pull requests have been filed to merge these equations into the main
**nephro** package
