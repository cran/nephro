\name{Convert.mgdL.mmolL}
\alias{Convert.mgdL.mmolL}
\title{
Convertion formula to change from mg/dL to mmol/L and vice versa
}
\description{
mg/dL <--> mmol/L conversion
}
\usage{
Convert.mgdL.mmolL(mg.dL, mmol.L)
## If your variable X is expressed in mg/dL, use: Convert.mgdL.mmolL(mg.dL=X)
## If your variable X is expressed in mmol/L, use: Convert.mgdL.mmolL(mmol.L=X)
}
\arguments{
  \item{mg.dL}{Provide a numeric vector in mg/dL to be converted into mmol/L}
  \item{mmol.L}{Provide a numeric vector in mmol/L to be converted into mg/dL}
}
\value{
A numeric vector with the converted values.
}
\author{
Cristian Pattaro
}

\seealso{
\code{\link[nephro]{Convert.mgdL.umolL}}
}
\keyword{ conversion }
\keyword{ creatinine }
