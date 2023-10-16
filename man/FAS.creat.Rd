\name{FAS.creat}
\alias{FAS.creat}
\title{
FAS equation for serum creatinine
}
\description{
Full age spectrum (FAS) equation to estimate GFR based on serum creatinine
}
\usage{
FAS.creat(creatinine, sex, age)
}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  }
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
   Pottel, H., \emph{et al}. An estimated glomerular filtration rate equation for the full age spectrum. \emph{Nephrol Dial Transplant.} 2016; \strong{5}: 798-806.}

\author{Janina Herold}
\seealso{
\code{\link{FAS.creat.cys}}, \code{\link{FAS.cys}}
}
\keyword{ FAS }
\keyword{ creatinine }
\keyword{ Pottel }
