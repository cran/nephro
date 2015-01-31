\name{Virga}
\alias{Virga}
\title{Virga's formula}
\description{
Virga's formula is based on serum creatinine, sex, age, and body weight. 
}
\usage{
Virga(creatinine, sex, age, wt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{creatinine}{Numeric vector with serum or plasma creatinine values in mg/dl}
  \item{sex}{Numeric 0/1 vector: 0 for females, 1 for males}
  \item{age}{Numeric vector with age in years}
  \item{wt}{Numeric vector with weight in kg}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
A numeric vector with eGFR values in \emph{ml/min/1.73 m^2} 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Virga G, \emph{et al}. A new equation for estimating renal function using age, body weight and serum creatinine. \emph{Nephron Clin Pract} 2007; \strong{105}: c43-53. 
}
\author{Cristian Pattaro}

%% \note{~~further notes~~}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[nephro]{CG}}, \code{\link{MDRD4}}
}

\keyword{ Virga }
\keyword{ creatinine }
\keyword{ clearance }
