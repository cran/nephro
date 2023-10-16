\name{EKFC.creat}
\alias{EKFC.creat}
\title{EKFC equation for serum creatinine}
\description{EKFC equation for serum creatinine modified from FAS equation
}
\usage{
EKFC.creat(creatinine, sex, age)
}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
}
\value{
The function returns a numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Pottel, H, \emph{et al}. Development and Validation of a Modified Full Age Spectrum Creatinine-Based Equation to Estimate Glomerular Filtration Rate : A Cross-sectional Analysis of Pooled Data. \emph{N Engl J Med} 2021; \strong{174}: 183-191.
}
\author{
Janina Herold
}

\seealso{
\code{\link{EKFC.cys}}, \code{\link{CKDEpi.creat}}
}
\keyword{ GFR }
\keyword{ EFKF }
\keyword{ creatinine }
