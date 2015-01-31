\name{CG}
\alias{CG}
\title{
Cockroft and Gault equation
}
\description{
Creatinine clearance is estimated with the Cockroft and Gault formula.
}
\usage{CG(creatinine, sex, age, wt)
}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{wt}{Numeric vector with weight in kg}
}
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
Cockroft DW, Gault MH. Prediction of creatinine clearance from serum creatinine. \emph{Nephron} 1976; \strong{16}: 31-41.\cr
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link[nephro]{CKDEpi.creat}}, \code{\link{MDRD4}}, \code{\link{Virga}}
}
\keyword{ MDRD }
\keyword{ creatinine }
\keyword{ clearance }
