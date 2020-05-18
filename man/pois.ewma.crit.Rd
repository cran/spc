\name{pois.ewma.crit}
\alias{pois.ewma.crit}
\title{Compute ARLs of Poisson EWMA control charts}
\description{Computation of the (zero-state) Average Run Length (ARL) at given mean \code{mu}.}
\usage{pois.ewma.crit(lambda, L0, mu0, z0, AU=3, sided="two", design="sym", rando=FALSE,
mcdesign="transfer", N=101, jmax=4)}
\arguments{
\item{lambda}{smoothing parameter of the EWMA p control chart.}
\item{L0}{value of the so-called in-control Average Run Length (ARL) for the Poisson EWMA control chart.}
\item{mu0}{in-control mean.}
\item{z0}{so-called headstart (give fast initial response).}
\item{AU}{in case of the lower chart deployed as reflecting upper barrier -- might be increased step by step until the resulting lower limit
does not change anymore.}
\item{sided}{distinguishes between one- and two-sided EWMA control chart by choosing \code{"upper"}, \code{"lower"}, and \code{"two"}, respectively.}
\item{design}{distinguishes between limits symmetric to the in-control mean \code{mu0} and an ARL-unbiased design (ARL maximum at \code{mu0}); use the shortcuts \code{"sym"} and \code{"unb"}, respectively, please.}
\item{rando}{Switch between the standard limit treatment, \code{FALSE},
and an additional randomisation (to allow `perfect' ARL calibration) by setting \code{TRUE}.
If randomisation is used, then the corresponding probailities, \code{gL} and \code{gU} are determined, appropriately.}
\item{mcdesign}{choose either \code{"classic"} which follows Borror, Champ and Rigdon (1998), or the more
sophisticated \code{"transfer"} which improves the accuracy heavily.}
\item{N}{number of states of the approximating Markov chain; is equal to the dimension of the resulting linear equation system.}
\item{jmax}{number of digits for the to be calculated factors \code{A} (sort of accuracy).}
}
\details{
The monitored data follow a Poisson distribution with \code{mu}.
Here we solve the inverse task to the usual ARL calculation. Hence, determine the control limit factors
so that the in-control ARL is (roughly) equal to \code{L0}.
The ARL values underneath the routine are determined by Markov chain approximation.
The algorithm is just a grid search that takes care of the discrete ARL behavior.
}
\value{Return one or two values being he control limit factors.}
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
## Borror, Champ and Rigdon (1998), page 30, original value is A = 2.8275
mu0 <- 4
lambda <- 0.2
L0 <- 351
A <- pois.ewma.crit(lambda, L0, mu0, mu0, mcdesign="classic")
print(round(A, digits=4))

## Morais and Knoth (2020), Table 2, lambda = 0.27 column
lambda <- 0.27
L0 <- 1233.4
ccgg <- pois.ewma.crit(lambda,1233.4,mu0,mu0,design="unb",rando=TRUE,mcdesign="transfer")
print(ccgg, digits=3)
}
\keyword{ts}
