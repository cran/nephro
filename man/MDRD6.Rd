\name{MDRD6}
\alias{MDRD6}
\title{
6-parameter MDRD study equation
}
\description{
GFR is estimated with the 6-parameter Modification of Diet in Renal Disease (MDRD) study equation.
}
\usage{
MDRD6(creatinine, sex, age, albumin, BUN, ethnicity, method = 'IDMS')
}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{albumin}{Numeric vector with serum or plasma albumin in g/dl}
  \item{BUN}{Numeric vector with blood urea nitrogen levels in mg/dl}
  \item{ethnicity}{Numeric vector with 0 for non-Black and 1 for Black individuals}
  \item{method}{Defaults is 'IDMS' for IDMS-traceable creatinine; write 'other' if not IDMS}
}
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
Levey AS, \emph{et al}. A more accurate method to estimate glomerular filtration rate from serum creatinine: a new prediction equation. Modification of Diet in Renal Disease Study Group. \emph{Ann Intern Med.} 1999;\strong{130}(6):461-70.\cr
Levey AS, \emph{et al}. Using standardized serum creatinine values in the modification of diet in renal disease study equation for estimating glomerular filtration rate. \emph{Ann Intern Med}. 2006;\strong{145}:247-254.
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link{MDRD4}}
}
\keyword{ ~MDRD }
\keyword{ ~creatinine }
\keyword{ ~albumin }
\keyword{ ~blood urea nitrogen }
