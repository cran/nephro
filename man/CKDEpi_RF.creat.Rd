\name{CKDEpi_RF.creat}
\alias{CKDEpi_RF.creat}
\title{Race-free CKD-EPI equation for serum creatinine
}
\description{
GFR is estimated with the CKD-EPI Study equation based on serum creatinine without the ethnicity coefficient.
}
\usage{CKDEpi_RF.creat(creatinine, sex, age)}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
}
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Inker LA, \emph{et al}. New creatinine- and cystatin C-based equations to estimate GFR without race. \emph{N Engl J Med} 2021; \strong{385}: 1737-1749.
}
\author{
Ryosuke Fujii
}

\seealso{
\code{\link{CKDEpi.creat}}, \code{\link{CKDEpi.creat.cys}}, \code{\link{CKDEpi.cys}}, \code{\link{CKDEpi_RF.creat.cys}}
}
\keyword{ CKD-EPI }
\keyword{ creatinine }
