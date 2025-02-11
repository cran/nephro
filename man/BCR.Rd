\name{BCR}
\alias{BCR}
\title{
BUN-to-creatinine ratio
}
\description{
Calculate the plasma or serum BUN-to-creatinine ratio
}
\usage{BCR(BUN, creatinine) 
}
\arguments{
  \item{BUN}{Numeric vector in mg/dL}
  \item{creatinine}{Numeric vector in mg/dL; 0s will be automatically converted into NAs}
}
\value{
A numeric vector with BCR values.
}
\references{
   Dossetor JB. Creatininemia Versus Uremia. The relative significance of blood urea nitrogen and serum creatinine concentrations in azotemia. \emph{Ann Int Med.} 1966; \strong{65}(6): 1287-1299.}
\author{
Cristian Pattaro
}

\seealso{
\code{\link[nephro]{UCR}}
}
\keyword{ BUN }
\keyword{ blood urea nitrogen }
\keyword{ urea }
\keyword{ creatinine }
\keyword{ BCR }
