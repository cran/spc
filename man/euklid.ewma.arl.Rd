\name{euklid.ewma.arl}
\alias{euklid.ewma.arl}
\title{Compute ARLs of Poisson NCS-EWMA control charts}
\description{Computation of the (zero-state) Average Run Length (ARL) at given Poisson mean \code{mu}.}
\usage{euklid.ewma.arl(gX, gY, kL, kU, mu, y0, r0=0)}
\arguments{
\item{gX}{first and}
\item{gY}{second integer forming the rational lambda = gX/(gX+gY), lambda mimics the usual EWMA smoothing constant.}
\item{kL}{lower control limit of the NCS-EWMA control chart, integer.}
\item{kU}{upper control limit of the NCS-EWMA control chart, integer.}
\item{mu}{mean value of Poisson distribution.}
\item{y0}{headstart like value -- it is proposed to use the in-control mean.}
\item{r0}{further element of the headstart -- deviating from the default should be done only in case of full understanding of the scheme.}
}
\details{
A new idea of applying EWMA smoothing to count data based on integer divison with remainders.
It is highly recommended to read the corresponding paper (see below).
}
\value{Return single value which resemble the ARL.}
\references{
A. C. Rakitzis, P. Castagliola, P. E. Maravelakis (2015),
A new memory-type monitoring technique for count data,
Computers and Industrial Engineering 85, 235-247.
}
\author{Sven Knoth}
\seealso{later.}
\examples{
# RCM (2015), Table 12, page 243, first NCS column
gX <- 5
gY <- 24
kL <- 16
kU <- 24
mu0 <- 20
#L0 <- euklid.ewma.arl(gX, gY, kL, kU, mu0, mu0)
# should be 1219.2
 
}
\keyword{ts}
