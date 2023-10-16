\name{Schwartz.Bedside}
\alias{Schwartz.Bedside}
\title{Bedside IDMS-traceable Schwartz GFR Calculator for Children}
\description{
GFR is estimated with the Bedside Schwartz equation for Children based on IDMS serum or plasma creatinine.\cr
This equation is valid in the 1-17 years age range.\cr
}
\usage{Schwartz.Bedside(creatinine, ht, age)}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{ht}{Numeric vector with height in cm}
  \item{age}{Numeric vector with age in years}
}
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Schwartz GJ, \emph{et al}. New equations to estimate GFR in children with CKD. \emph{J Am Soc Nephrol} 2009; \strong{20}: 629-637.
}
\author{
Andrew Srisuwananukorn 
}

\keyword{ Schwartz }
\keyword{ creatinine }
