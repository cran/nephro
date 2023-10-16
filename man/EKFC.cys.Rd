\name{EKFC.cys}
\alias{EKFC.cys}
\title{EKFC equation for cystatin C}
\description{EKFC equation for cystatin C that includes the sex coefficient, as proposed by Pottel et al., \emph{N Engl J Med} 2023
}
\usage{
EKFC.cys(cystatin, sex, age)
}
\arguments{
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/l}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
}
\value{
The function returns a numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Pottel, H, \emph{et al}. Cystatin C-Based Equation to Estimate GFR without the Inclusion of Race and Sex. \emph{N Engl J Med} 2023; \strong{388}: 333-343.
}
\author{
Janina Herold
}

\seealso{
\code{\link{EKFC_SF.cys}}
}
\keyword{ GFR }
\keyword{ EFKF }
\keyword{ cystatin }
\keyword{ sex }
