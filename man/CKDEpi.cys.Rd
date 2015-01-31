\name{CKDEpi.cys}
\alias{CKDEpi.cys}
\title{CKD-EPI equation for cystatin C}
\description{
GFR is estimated with the CKD-EPI equation for cystatin C proposed by Inker et al., N Engl J Med 2012
}
\usage{
CKDEpi.cys(cystatin, sex, age)
}
\arguments{
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
}
\value{
The function returns a numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Inker LA, \emph{et al}. Estimating glomerular filtration rate from serum creatinine and cystatin C. \emph{N Engl J Med} 2012; \strong{367}: 20-29.
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link{CKDEpi.creat}}, \code{\link{CKDEpi.creat.cys}}
}
\keyword{ GFR }
\keyword{ CKD-EPI }
\keyword{ creatinine }
\keyword{ cystatin }
