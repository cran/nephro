\name{UPCR}
\alias{UPCR}
\title{
Urinary protein-to-creatinine ratio
}
\description{
Calculate the urinary protein-to-creatinine ratio
}
\usage{UPCR(u.protein, u.creatinine) 
}
\arguments{
  \item{u.protein}{Numeric vector in mg/dL}
  \item{u.creatinine}{Numeric vector in mg/dL; 0s will be automatically converted into NAs}
}
\value{
A numeric vector with UPCR values.
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link[nephro]{UACR}}
}
\keyword{ protein }
\keyword{ creatinine }
\keyword{ UPCR }
\keyword{ proteinuria }
