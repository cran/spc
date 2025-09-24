\name{pois.ewma.arl}
\alias{pois.ewma.arl}
\title{Compute ARLs of Poisson EWMA control charts}
\description{Computation of the (zero-state) Average Run Length (ARL) at given mean \code{mu}.}
\usage{pois.ewma.arl(lambda, AL, AU, mu0, z0, mu, sided="two", rando=FALSE, gL=0, gU=0,
mcdesign="transfer", N=101, OLD=FALSE)}
\arguments{
\item{lambda}{smoothing parameter of the EWMA p control chart.}
\item{AL, AU}{factors to build the lower and upper control limit, respectively, of the Poisson EWMA control chart.}
\item{mu0}{in-control mean.}
\item{z0}{so-called headstart (give fast initial response).}
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
\item{OLD}{Switch between the original implemenation, used in Morais and Knoth (2020), and a slightly modified new one (more accurate).}
}
\details{
The monitored data follow a Poisson distribution with \code{mu}.
The ARL values of the resulting EWMA control chart are determined by Markov chain approximation.
We follow the algorithm given in Borror, Champ and Rigdon (1998).
However, by setting \code{mcdesign="transfer"} (now the default) from Morais and Knoth (2020),
the accuracy is considerably improved.
}
\value{Return single value which resembles the ARL.}
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
ARL1 <- rep(NA, length(mu1))
for ( i in 1:length(mu1) )
  ARL1[i] <- pois.ewma.arl(lambda, A, A, mu0, mu0, mu1[i], mcdesign="classic")
print(cbind(mu1, round(ARL1, digits=1)))

## the same numbers with improved accuracy
ARL2 <- rep(NA, length(mu1))
for ( i in 1:length(mu1) )
  ARL2[i] <- pois.ewma.arl(lambda, A, A, mu0, mu0, mu1[i], mcdesign="transfer")
print(cbind(mu1, round(ARL2, digits=1)))

## Morais and Knoth (2020), Table 2, lambda = 0.27 column
lambda <- 0.27
AL <- 3.0870
AU <- 3.4870
gL <- 0.001029
gU <- 0.000765
mu0 <- 20
mu1  <- c(16, 18, 19.99, mu0, 20.01, 22, 24)
ARL3 <- rep(NA, length(mu1))
for ( i in 1:length(mu1) )
  ARL3[i] <- pois.ewma.arl(lambda,AL,AU,mu0,mu0,mu1[i],rando=TRUE,gL=gL,gU=gU, N=101, OLD=TRUE)
print(cbind(mu1, round(ARL3, digits=1)))
}
\keyword{ts}
