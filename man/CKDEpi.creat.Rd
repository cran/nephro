\name{CKDEpi.creat}
\alias{CKDEpi.creat}
\title{CKD-EPI equation for serum creatinine}
\description{
GFR is estimated with the CKD-EPI Study equation based on IDMS serum or plasma creatinine.\cr
}
\usage{CKDEpi.creat(creatinine, sex, age, ethnicity)}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{ethnicity}{Numeric vector with 0 for non-Black and 1 for Black individuals}
}
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Inker LA, \emph{et al}. Estimating glomerular filtration rate from serum creatinine and cystatin C. \emph{N Engl J Med} 2012; \strong{367}: 20-29.
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link{CKDEpi.creat.cys}}, \code{\link{CKDEpi.cys}}, \code{\link{CKDEpi_RF.creat}}
}
\keyword{ CKD-EPI }
\keyword{ creatinine }
