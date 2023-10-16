\name{EKFC_SF.cys}
\alias{EKFC_SF.cys}
\title{Sex-free EKFC equation for cystatin C}
\description{EKFC equation for cystatin C without the sex coefficient as proposed by Pottel et al., \emph{N Engl J Med} 2023
}
\usage{
EKFC_SF.cys(cystatin, age)
}
\arguments{
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/l}
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
\code{\link{EKFC.cys}}
}
\keyword{ GFR }
\keyword{ EFKF }
\keyword{ cystatin }
