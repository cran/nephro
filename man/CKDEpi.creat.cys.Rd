\name{CKDEpi.creat.cys}
\alias{CKDEpi.creat.cys}
\title{
CKD-EPI equation for creatinine and cystatin C
}
\description{
CKD-EPI equation to estimate GFR based on a combination of creatinine and cystatin C
}
\usage{
CKDEpi.creat.cys(creatinine, cystatin, sex, age, ethnicity)
}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{ethnicity}{Numeric vector with 0 for non-Black and 1 for Black individuals}
}
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Inker LA, \emph{et al}. Estimating glomerular filtration rate from serum creatinine and cystatin C. \emph{N Engl J Med} 2012; \strong{367}: 20-29.}
\author{Cristian Pattaro}

\seealso{
\code{\link{CKDEpi.creat}}, \code{\link{CKDEpi.cys}}
}
\keyword{ CKD-EPI }
\keyword{ creatinine }
\keyword{ cystatin }
