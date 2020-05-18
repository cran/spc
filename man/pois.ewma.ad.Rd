\name{pois.ewma.ad}
\alias{pois.ewma.ad}
\title{Compute steady-state ARLs of Poisson EWMA control charts}
\description{Computation of the steady-state Average Run Length (ARL) at given mean \code{mu}.}
\usage{pois.ewma.ad(lambda, AL, AU, mu0, mu, sided="two", rando=FALSE, gL=0, gU=0,
mcdesign="classic", N=101)}
\arguments{
\item{lambda}{smoothing parameter of the EWMA p control chart.}
\item{AL, AU}{factors to build the lower and upper control limit, respectively, of the Poisson EWMA control chart.}
\item{mu0}{in-control mean.}
\item{mu}{actual mean.}
\item{sided}{distinguishes between one- and two-sided EWMA control chart by choosing
\code{"upper"}, \code{"lower"}, and \code{"two"}, and \code{"zwei"}, respectively.}
\item{rando}{Switch between the standard limit treatment, \code{FALSE},
and an additional randomisation (to allow `perfect' ARL calibration) by setting \code{TRUE}.
If randomisation is used, then set the corresponding probailities, \code{gL} and \code{gU}, appropriately.}
\item{gL, gU}{If the EWMA statistic is at the limit (approximately), then an alarm is triggered with probability
\code{gL} and \code{gU} for the lower and upper limit, respectively.}
\item{mcdesign}{choose either \code{"classic"} which follows Borror, Champ and Rigdon (1998), or the more
sophisticated \code{"transfer"} which improves the accuracy heavily.}
\item{N}{number of states of the approximating Markov chain; is equal to the dimension of the
resulting linear equation system.}
}
\details{
The monitored data follow a Poisson distribution with \code{mu}.
The ARL values of the resulting EWMA control chart are determined by Markov chain approximation.
We follow the algorithm given in Borror, Champ and Rigdon (1998). The function is in an early development phase.
}
\value{Return single value which resembles the steady-state ARL.}
\references{
C. M. Borror, C. W. Champ and S. E. Rigdon (1998)
Poisson EWMA control charts,
\emph{Journal of Quality Technonlogy} 30(4), 352-361.

M. C. Morais and S. Knoth (2020)
Improving the ARL profile and the accuracy of its calculation for Poisson EWMA charts,
\emph{Quality and Reliability Engineering International} 36(3), 876-889.
}
\author{Sven Knoth}
\seealso{later.}
\examples{
## Borror, Champ and Rigdon (1998), Table 2, PEWMA column
mu0 <- 20
lambda <- 0.27
A <- 3.319
mu1  <- c(2*(3:15), 35)
ARL1 <- AD1 <- rep(NA, length(mu1))
for ( i in 1:length(mu1) ) {
  ARL1[i] <- round(pois.ewma.arl(lambda,A,A,mu0,mu0,mu1[i],mcdesign="classic"),digits=1)
  AD1[i]  <- round(pois.ewma.ad(lambda,A,A,mu0,mu1[i],mcdesign="classic"),digits=1)
}
print( cbind(mu1, ARL1, AD1) )

## Morais and Knoth (2020), Table 2, lambda = 0.27 column
## randomisation not implemented for pois.ewma.ad()
lambda <- 0.27
AL <- 3.0870
AU <- 3.4870
gL <- 0.001029
gU <- 0.000765
mu2  <- c(16, 18, 19.99, mu0, 20.01, 22, 24)
ARL2 <- AD2 <- rep(NA, length(mu2))
for ( i in 1:length(mu2) ) {
  ARL2[i] <- round(pois.ewma.arl(lambda,AL,AU,mu0,mu0,mu2[i],rando=FALSE), digits=1)
  AD2[i] <- round(pois.ewma.ad(lambda,AL,AU,mu0,mu2[i],rando=FALSE), digits=1)
}
print( cbind(mu2, ARL2, AD2) )
}
\keyword{ts}
