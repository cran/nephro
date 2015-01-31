\name{Stevens.cys2}
\alias{Stevens.cys2}
\title{Stevens' formula for serum cystatin C, age, and sex}
\description{
GFR is estimated with the 2nd formula proposed by Stevens et al. (Am J Kidney Dis 2008), where cystatin C is weighted by sex and age
}
\usage{Stevens.cys2(cystatin, sex, age, ethnicity)}
\arguments{
  \item{cystatin}{Numeric vector with serum or plasma cystatin C values in mg/dl}
  \item{sex}{Numeric vector with 0 for females and 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{ethnicity}{Numeric vector with 0 for non-Black and 1 for Black individuals}
}
\value{A numeric vector with eGFR values in ml/min/1.73 \eqn{m^2}.}
\references{
Stevens LA, \emph{et al}. Estimating GFR using serum cystatin C alone and in combination with serum creatinine: a pooled analysis of 3,418 individuals with CKD. \emph{Am J Kidney Dis} 2008; \strong{51}: 395-406. 
}
\author{Cristian Pattaro}
\seealso{
\code{\link{Stevens.cys1}}, \code{\link{Stevens.creat.cys}}, \code{\link{CKDEpi.cys}}
}
\keyword{ GFR }
\keyword{ cystatin }
