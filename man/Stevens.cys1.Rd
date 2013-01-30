\name{Stevens.cys1}
\alias{Stevens.cys1}
\title{
GFR estimation using serum cystatin C
}
\description{
GFR is estimated with the 1st formula proposed by Stevens et al. (Am J Kidney Dis 2008), i.e.: as a simple transformation of cystatin C, without using any other information
}
\usage{
Stevens.cys1(cystatin)
}
\arguments{
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/dl}
}
\value{
A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
Stevens LA, \emph{et al}. Estimating GFR using serum cystatin C alone and in combination with serum creatinine: a pooled analysis of 3,418 individuals with CKD. \emph{Am J Kidney Dis} 2008;\strong{51}:395-406. 
}
\author{
Cristian Pattaro
}
\seealso{
\code{\link{Stevens.cys2}, \link{Stevens.creat.cys}, \link{CKDEpi.cys}}
}
\keyword{ ~GFR }
\keyword{ ~cystatin }
