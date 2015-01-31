\name{nephro-package}
\alias{nephro-package}
\alias{nephro}
\docType{package}
\title{
Biostatistics Utilities for Nephrology
}
\description{
This package contains a set of tools for the estimation of renal function. Renal function is assessed by means of the Glomerular Filtration Rate (GFR) which can be estimated using different biomarkers. The most commonly used ones are serum or plasma creatinine and cystatin C. Included in this package are the following GFR estimating functions: the Modification of Diet in Renal Disease (MDRD) study equation based on four (\code{\link[nephro]{MDRD4}}) or six (\code{\link[nephro]{MDRD6}}) parameters (see Levey 1999, 2006); the CKD-Epi equations for creatinine (\code{\link[nephro]{CKDEpi.creat}}), cystatin C (\code{\link[nephro]{CKDEpi.cys}}), and a combination of creatinine and cystatin C (\code{\link[nephro]{CKDEpi.creat.cys}}) (see Inker 2012); the three equations proposed by Stevens 2008 based on cystatin C only (\code{\link[nephro]{Stevens.cys1}}), age- and sex-weighted cystatin C (\code{\link[nephro]{Stevens.cys2}}), and a combination of cystatin C and creatinine (\code{\link[nephro]{Stevens.creat.cys}}); the classic Cockroft and Gault 1976 equation for creatinine clearance estimation (\code{\link[nephro]{CG}}); and the recent equation by Virga (2007) (\code{\link[nephro]{Virga}}). A comparative description of several functions included in the package can be found in Pattaro (2013). An extensive literature does exist that compares the methods described.
}
\details{
\tabular{ll}{
Package: \tab nephro\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-01-31\cr
License: \tab GPLv3\cr
}
}
\author{
Cristian Pattaro <cristian.pattaro@eurac.edu> \cr

Acknowledgements: Max Plischke (contributed examples and bug reporting)
}

\references{
  Citing this package:\cr
  - Pattaro C, Riegler P, Stifter G, Modenese M, Minelli C, Pramstaller PP. Estimating the glomerular filtration rate in the general population using different equations: effects on classification and association. \emph{Nephron Clin Pract} 2013; \strong{123}(1-2):102-11.\cr
  
  Formulas:\cr
  - Cockroft DW, Gault MH. Prediction of creatinine clearance from serum creatinine. \emph{Nephron} 1976; \strong{16}: 31-41.\cr
  - Inker LA, \emph{et al}. Estimating glomerular filtration rate from serum creatinine and cystatin C. \emph{N Engl J Med} 2012; \strong{367}: 20-9. \cr
  - Levey AS, \emph{et al}. A more accurate method to estimate glomerular filtration rate from serum creatinine: a new prediction equation. Modification of Diet in Renal Disease Study Group. \emph{Ann Intern Med}. 1999; \strong{130}(6): 461-70.\cr
  - Levey AS, \emph{et al}. Using standardized serum creatinine values in the modification of diet in renal disease study equation for estimating glomerular filtration rate. \emph{Ann Intern Med.} 2006; \strong{145}: 247-54. \cr
  - Stevens LA, \emph{et al}. Estimating GFR using serum cystatin C alone and in combination with serum creatinine: a pooled analysis of 3,418 individuals with CKD. \emph{Am J Kidney Dis} 2008; \strong{51}: 395-406. \cr
  - Virga G, \emph{et al}. A new equation for estimating renal function using age, body weight and serum creatinine. \emph{Nephron Clin Pract} 2007; \strong{105}: c43-53.\cr 
  
  On IDMS calibration:\cr
  - Levey AS, \emph{et al}. Expressing the Modification of Diet in Renal Disease Study equation for estimating glomerular filtration rate with standardized serum creatinine values. \emph{Clin Chem} 2007; \strong{53}:766-72.\cr
  - Matsushita K, \emph{et al}. Comparison of risk prediction using the CKD-EPI equation and the MDRD study equation for estimated glomerular filtration rate. \emph{J Am Med Assoc} 2012; \strong{307}:1941-51.\cr
  - Skali H, \emph{et al}. Prognostic assessment of estimated glomerular filtration rate by the new Chronic Kidney Disease Epidemiology Collaboration equation in comparison with the Modification of Diet in Renal Disease Study equation. \emph{Am Heart J} 2011; \strong{162}:548-54.\cr

}

\keyword{
GFR, renal function, kidney, MDRD, CKD-EPI, creatinine, cystatin
}

\examples{
# Comparison between different equations

creat <- c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
cyst  <- c(1.1, 0.95, 1.1, 1.0, 1.3, 1.2)
sex <- c(1, 1, 1, 0, 0, 0)
age <- c(60, 65, 43, 82, 71, 55)
ethn <- round(runif(6))
wt <- c(70, 80, 60, 55, 87, 71)

eGFR <- data.frame(creat, cyst)
eGFR$MDRD4 <- MDRD4(creat, sex, age, ethn, 'IDMS') 
eGFR$CKDEpi.creat <- CKDEpi.creat(creat, sex, age, ethn) 
eGFR$CKDEpi.cys <- CKDEpi.cys(cyst, sex, age)
eGFR$CKDEpi.creat.cys <- CKDEpi.creat.cys(creat, cyst, sex, age, ethn)
eGFR$Stevens.cys1 <- Stevens.cys1(cyst)
eGFR$Stevens.cys2 <- Stevens.cys2(cyst, sex, age, ethn)
eGFR$Stevens.creat.cys <- Stevens.creat.cys(creat, cyst, sex, age, ethn)
eGFR$cg <- CG(creat, sex, age, wt)
eGFR$virga <- Virga(creat, sex, age, wt)

pairs(eGFR[,3:11])


# For use with non-IDMS calibrated creatinine 
# several authors (see references) suggested
# a 5\% creatinine adjustment

creat <- c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
sex <- c(1, 1, 1, 0, 0, 0)
age <- c(60, 65, 43, 82, 71, 55)
ethn <- round(runif(6))
gfr <- CKDEpi.creat(0.95*creat, sex, age, ethn) 
}
