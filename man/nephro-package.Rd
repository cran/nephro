\name{nephro-package}
\alias{nephro-package}
\alias{nephro}
\docType{package}
\title{
Biostatistics utilities for nephrology
}
\description{
Collection of functions to estimate the Glomerular Filtration Rate (GFR) using different biomarkers and other phenotypes of interest for nephrology.
}
\details{
\tabular{ll}{
Package: \tab nephro\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-30\cr
License: \tab GPLv3\cr
}
}
\author{
Cristian Pattaro <cristian.pattaro@eurac.edu>
}
\references{
   Inker LA, \emph{et al}. Estimating glomerular filtration rate from serum creatinine and cystatin C. \emph{N Engl J Med} 2012;\strong{367}:20-29. \cr 

   Levey AS, \emph{et al}. A more accurate method to estimate glomerular filtration rate from serum creatinine: a new prediction equation. Modification of Diet in Renal Disease Study Group. \emph{Ann Intern Med} 1999;\strong{130}(6):461-470. \cr

   Levey AS, \emph{et al}. Using standardized serum creatinine values in the modification of diet in renal disease study equation for estimating glomerular filtration rate. \emph{Ann Intern Med} 2006;\strong{145}:247-254. \cr

   Pattaro C, \emph{et al}. Estimating the glomerular filtration rate in the general population using different equations: effects on classification and association. \emph{Submitted}

   Stevens LA, \emph{et al}. Estimating GFR using serum cystatin C alone and in combination with serum creatinine: a pooled analysis of 3,418 individuals with CKD. \emph{Am J Kidney Dis} 2008;\strong{51}:395-406. \cr

   Virga G, \emph{et al}. A new equation for estimating renal function using age, body weight and serum creatinine. \emph{Nephron Clin Pract} 2007;\strong{105}:c43-53. \cr
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

eGFR <- data.frame(creat, cyst)
eGFR$MDRD4 <- MDRD4(creat, sex, age, ethn, 'IDMS') 
eGFR$CKDEpi.creat <- CKDEpi.creat(creat, sex, age, ethn) 
eGFR$CKDEpi.cys <- CKDEpi.cys(cyst, sex, age)
eGFR$CKDEpi.creat.cys <- CKDEpi.creat.cys(creat, cyst, sex, age, ethn)
eGFR$Stevens.cys1 <- Stevens.cys1(cyst)
eGFR$Stevens.cys2 <- Stevens.cys2(cyst, sex, age, ethn)
eGFR$Stevens.creat.cys <- Stevens.creat.cys(creat, cyst, sex, age, ethn)

pairs(eGFR[,3:9])

}
