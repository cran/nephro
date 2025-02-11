\name{nephro-package}
\alias{nephro-package}
\alias{nephro}
\docType{package}
\title{
Biostatistics utilities for nephrology
}
\description{
Set of tools for the estimation of kidney function and related traits. 

Included are the following glomerular filtration rate (GFR) estimating functions: 
the Modification of Diet in Renal Disease (MDRD) study equations based on four (\code{\link[nephro]{MDRD4}}) or six (\code{\link[nephro]{MDRD6}}) parameters (Levey 1999; Levey 2006);
the CKD-Epi equations for serum creatinine with the race coefficient (\code{\link[nephro]{CKDEpi.creat}}) and without the race coefficient (\code{\link[nephro]{CKDEpi2021.creat}});
the CKD-Epi equation for cystatin C (\code{\link[nephro]{CKDEpi.cys}}); the CKD-Epi equation for the combination of creatinine and cystatin C with (\code{\link[nephro]{CKDEpi.creat.cys}}) and without (\code{\link[nephro]{CKDEpi2021.creat.cys}}) the race coefficient (Inker 2012; Inker 2021);
the three equations proposed by Stevens 2008 based on cystatin C only (\code{\link[nephro]{Stevens.cys1}}), age- and sex-weighted cystatin C (\code{\link[nephro]{Stevens.cys2}}), and a combination of cystatin C and creatinine (\code{\link[nephro]{Stevens.creat.cys}});
the classic Cockroft and Gault 1976 equation for creatinine clearance estimation (\code{\link[nephro]{CG}}); the equation by Virga (2007) (\code{\link[nephro]{Virga}});
the race-free equations developed by the European Kidney Function Consortium (EKFC) including sex and age based on serum creatinine (\code{\link[nephro]{EKFC.creat}}) (Pottel 2021) and based on serum cystatin C with (\code{\link[nephro]{EKFC.cys}}) and without (\code{\link[nephro]{EKFC_SF.cys}}) the sex coefficient (Pottel 2023);
the full age spectrum (FAS) equations using serum creatinine (\code{\link[nephro]{FAS.creat}}) (Pottel 2016), cystatin C (\code{\link[nephro]{FAS.cys}}), and their combination (\code{\link[nephro]{FAS.creat.cys}}) (Pottel 2017);
the Schwartz bedside formula (\code{\link[nephro]{Schwartz.Bedside}}) (Schwartz 2009).

A comparative description of several GFR-estimating functions included in the initial version of the package can be found in Pattaro (2013). Extensive literature does exist that compares the methods described.

Additional serum biomarkers are the blood urea nitrogen (BUN)-to-creatinine ratio (\code{\link[nephro]{BCR}} and the urea-to-creatinine ratio \code{\link[nephro]{UCR}}.

Urine biomarkers included are the albumin-to-creatinine ratio \code{\link[nephro]{UACR}} and the protein-to-creatinine ratio \code{\link[nephro]{UPCR}}

Some unit conversion functions are provided to facilitate operations.
}

\details{
\tabular{ll}{
Package: \tab nephroext\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2025-30-01\cr
License: \tab GPLv3\cr
}
}
\author{
Cristian Pattaro <cristian.pattaro@eurac.edu>, Ryosuke Fujii <ryosuke.fujii@eurac.edu>, Janina Herold <janina.herold@klinik.uni-regensburg.de>\cr

Acknowledgements: Max Plischke, Joao Sabino, Xiao Jiang for bug reporting; Alexander Tolios, Thomas Winkler for input on the software package; Andrew Srisuwananukorn for contributing scripts; the 'nephro' user community for continuous feedback.
}

\references{
  Citing this package:\cr
  - Pattaro C, Riegler P, Stifter G, Modenese M, Minelli C, Pramstaller PP. Estimating the glomerular filtration rate in the general population using different equations: effects on classification and association. \emph{Nephron Clin Pract} 2013; \strong{123}(1-2):102-11.\cr
  
  Formulas:\cr
  - Cockroft DW, Gault MH. Prediction of creatinine clearance from serum creatinine. \emph{Nephron} 1976; \strong{16}: 31-41.\cr
  - Inker LA, \emph{et al}. Estimating glomerular filtration rate from serum creatinine and cystatin C. \emph{N Engl J Med} 2012; \strong{367}: 20-9. \cr
  - Inker LA, \emph{et al}. New Creatinine- and Cystatin C-based Equations to Estimate GFR without Race. \emph{N Engl J Med} 2021; \strong{385}: 1737-1749.\cr
  - Levey AS, \emph{et al}. A more accurate method to estimate glomerular filtration rate from serum creatinine: a new prediction equation. Modification of Diet in Renal Disease Study Group. \emph{Ann Intern Med}. 1999; \strong{130}(6): 461-70.\cr
  - Levey AS, \emph{et al}. Using standardized serum creatinine values in the modification of diet in renal disease study equation for estimating glomerular filtration rate. \emph{Ann Intern Med.} 2006; \strong{145}: 247-54. \cr
  - Stevens LA, \emph{et al}. Estimating GFR using serum cystatin C alone and in combination with serum creatinine: a pooled analysis of 3,418 individuals with CKD. \emph{Am J Kidney Dis} 2008; \strong{51}: 395-406. \cr
  - Virga G, \emph{et al}. A new equation for estimating renal function using age, body weight and serum creatinine. \emph{Nephron Clin Pract} 2007; \strong{105}: c43-53.\cr 
  - Pottel H, \emph{et al}. Development and Validation of a Modified Full Age Spectrum Creatinine-Based Equation to Estimate Glomerular Filtration Rate: A Cross-sectional Analysis of Pooled Data. \emph{Ann Intern Med} 2021; \strong{174}: 183-191.\cr
  - Pottel H, \emph{et al}. Cystatin C-Based Equation to Estimate GFR without the Inclusion of Race and Sex. \emph{N Engl J Med} 2023; \strong{388}: 333-343.\cr
  - Pottel H, \emph{et al}. An estimated glomerular filtration rate equation for the full age spectrum. \emph{Nephrol Dial Transplant} 2016; \strong{31}:798-806.\cr
  - Pottel H, \emph{et al}. Estimating glomerular filtration rate for the full age spectrum from serum creatinine and cystatin C \emph{Nephrol Dial Transplant} 2017; \strong{32}: 497-507.\cr
  - Schwartz GJ, \emph{et al}. New equations to estimate GFR in children with CKD. \emph{J Am Soc Nephrol} 2009; \strong{20}:629-637.\cr

  On IDMS calibration:\cr
  - Levey AS, \emph{et al}. Expressing the Modification of Diet in Renal Disease Study equation for estimating glomerular filtration rate with standardized serum creatinine values. \emph{Clin Chem} 2007; \strong{53}:766-72.\cr
  - Matsushita K, \emph{et al}. Comparison of risk prediction using the CKD-EPI equation and the MDRD study equation for estimated glomerular filtration rate. \emph{J Am Med Assoc} 2012; \strong{307}:1941-51.\cr
  - Skali H, \emph{et al}. Prognostic assessment of estimated glomerular filtration rate by the new Chronic Kidney Disease Epidemiology Collaboration equation in comparison with the Modification of Diet in Renal Disease Study equation. \emph{Am Heart J} 2011; \strong{162}:548-54.\cr

}

\concept{kidney function}
\concept{GFR}
\concept{CKD-EPI}
\concept{FAS}
\concept{EKFC}
\concept{MDRD}
\concept{creatinine}
\concept{cystatin}
\concept{UACR}
\concept{albuminuria}
\concept{proteinuria}
\concept{azotaemia}

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
eGFR$CKDEpi2021.creat <- CKDEpi2021.creat(creat, sex, age)
eGFR$CKDEpi.cys <- CKDEpi.cys(cyst, sex, age)
eGFR$CKDEpi.creat.cys <- CKDEpi.creat.cys(creat, cyst, sex, age, ethn)
eGFR$CKDEpi2021.creat.cys <- CKDEpi2021.creat.cys(creat, cyst, sex, age)
eGFR$Stevens.cys1 <- Stevens.cys1(cyst)
eGFR$Stevens.cys2 <- Stevens.cys2(cyst, sex, age, ethn)
eGFR$Stevens.creat.cys <- Stevens.creat.cys(creat, cyst, sex, age, ethn)
eGFR$cg <- CG(creat, sex, age, wt)
eGFR$virga <- Virga(creat, sex, age, wt)

pairs(eGFR[,3:13])

# For use with non-IDMS calibrated creatinine 
# several authors (see references) suggested
# a 5\% creatinine adjustment

creat <- c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
sex <- c(1, 1, 1, 0, 0, 0)
age <- c(60, 65, 43, 82, 71, 55)
ethn <- round(runif(6))
gfr <- CKDEpi.creat(0.95*creat, sex, age, ethn) 
}