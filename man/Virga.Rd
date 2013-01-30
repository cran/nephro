\name{Virga}
\alias{Virga}
\title{
Virga's formula
}
\description{
Virga's formula is based on serum creatinine, sex, age, and body weight.
}
\usage{
Virga(creatinine, sex, age, wt)
}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{wt}{Numeric vector with weight in kg}
}
\value{
A numeric vector with eGFR values in \eqn{ml/min/1.73 m^2}.
}

\references{
   Virga G, \emph{et al}. A new equation for estimating renal function using age, body weight and serum creatinine. \emph{Nephron Clin Pract} 2007;\strong{105}:c43-53.
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link{CKDEpi.creat}}, \code{\link{MDRD4}}
}
\keyword{ ~GFR }
\keyword{ ~creatinine }
