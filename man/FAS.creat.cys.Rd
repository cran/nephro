\name{FAS.creat.cys}
\alias{FAS.creat.cys}
\title{FAS equation for creatinine and cystatin C}
\description{Full age spectrum (FAS) GFR estimation based on serum creatinine and cystatin C
}
\usage{
FAS.creat.cys(creatinine, cystatin, sex,age)
}

\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/l}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
}
\value{
The function returns a numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Pottel, H., \emph{et al}. An estimated glomerular filtration rate equation for the full age spectrum from serum creatinine and cystatin C. \emph{Nephrol Dial Transplant.} 2017; \strong{32}: 497-507.}
   
\author{Janina Herold}

\seealso{
\code{\link{FAS.cys}}
}
\keyword{ GFR }
\keyword{ EFKF }
\keyword{ cystatin }
