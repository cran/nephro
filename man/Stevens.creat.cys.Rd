\name{Stevens.creat.cys}
\alias{Stevens.creat.cys}
\title{Stevens' formula for a combination of serum creatinine and cystatin C}
\description{
GFR estimation using the 3rd formula proposed by Stevens et al. (Am J Kidney Dis 2008), which combines creatinine and cystatin C
}
\usage{
Stevens.creat.cys(creatinine, cystatin, sex, age, ethnicity)
}
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{ethnicity}{Numeric vector with 0 for non-Black and 1 for Black individuals}
}
\value{
The function returns a numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.
}
\references{
Stevens LA, \emph{et al}. Estimating GFR using serum cystatin C alone and in combination with serum creatinine: a pooled analysis of 3,418 individuals with CKD. \emph{Am J Kidney Dis} 2008;\strong{51}:395-406. 
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link{CKDEpi.creat.cys}}
}
\keyword{ ~GFR }
\keyword{ ~creatinine }
\keyword{ ~cystatin }
