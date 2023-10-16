\name{FAS.cys}
\alias{FAS.cys}
\title{
FAS equation for cystatin C
}
\description{
Full age spectrum (FAS) GFR estimation based on cystatin C
}
\usage{
FAS.cys(cystatin, sex, age)
}
\arguments{
  \item{cystatin}{Numeric vector with serum or plasma cystatin values in mg/l}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  }
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Pottel, H., \emph{et al}. An estimated glomerular filtration rate equation for the full age spectrum from serum creatinine and cystatin C. \emph{Nephrol Dial Transplant.} 2017; \strong{32}: 497-507.}

\author{Janina Herold}
\seealso{
\code{\link{FAS.creat}}, \code{\link{FAS.creat.cys}}
}
\keyword{ FAS }
\keyword{ cystatin }
\keyword{ Pottel }
