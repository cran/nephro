\name{UACR}
\alias{UACR}
\title{
Urinary albumin-to-creatinine ratio
}
\description{
Calculate the urinary albumin-to-creatinine ratio
}
\usage{UACR(u.albumin, u.creatinine) 
}
\arguments{
  \item{u.albumin}{Numeric vector in mg/dL}
  \item{u.creatinine}{Numeric vector in mg/dL; 0s will be automatically converted into NAs}
}
\value{
A numeric vector with UACR values.
}
\references{
Barratt TM, McLaine PN, Soothill JF. Albumin excretion as a measure of glomerular dysfunction in children. \emph{Arch Dis Child.} 1970. \strong{45}(242): 496-501
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link[nephro]{UPCR}}
}
\keyword{ albumin }
\keyword{ creatinine }
\keyword{ UACR }
\keyword{ albuminuria }
