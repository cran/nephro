\name{UCR}
\alias{UCR}
\title{
Urea-to-creatinine ratio
}
\description{
Calculate the plasma or serum urea-to-creatinine ratio
}
\usage{UCR(urea, creatinine) 
}
\arguments{
  \item{urea}{Numeric vector in mmol/L}
  \item{creatinine}{Numeric vector in umol/L; 0s will be automatically converted into NAs}
}
\value{
A numeric vector with UCR values.
}
\references{
   Kahn S, \emph{et al}. The significance of serum creatinine and the blood urea/serum creatinine ratio in azotaemia \emph{S. Afr. Med. J.} 1972; \strong{46}: 1828-1832.}
\author{
Cristian Pattaro
}

\seealso{
\code{\link[nephro]{BCR}}
}
\keyword{ urea }
\keyword{ creatinine }
\keyword{ UCR }
