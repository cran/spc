\name{xtewma.arl}
\alias{xtewma.arl}
\title{Compute ARLs of EWMA control charts, t distributed data}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of EWMA control charts monitoring the mean of t distributed data.}
\usage{xtewma.arl(l,c,df,mu,zr=0,hs=0,sided="two",limits="fix",mode="tan",q=1,r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{df}{degrees of freedom -- parameter of the t distribution.}
\item{mu}{true mean.}
\item{zr}{reflection border for the one-sided chart.}
\item{hs}{so-called headstart (enables fast initial response).}
\item{sided}{distinguishes between one- and two-sided EWMA control chart 
by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguishes between different control limits behavior.}
\item{mode}{Controls the type of variables substitution that might improve the numerical performance. Currently,
\code{"identity"}, \code{"sin"}, \code{"sinh"}, and \code{"tan"} (default) are provided.}
\item{q}{change point position. For \eqn{q=1} and
\eqn{\mu=\mu_0} and \eqn{\mu=\mu_1}, the usual
zero-state ARLs for the in-control and out-of-control case, respectively,
are calculated. For \eqn{q>1} and \eqn{\mu!=0} conditional delays, that is,
\eqn{E_q(L-q+1|L\ge q)}, will be determined.
Note that mu0=0 is implicitely fixed.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
}
\details{
In case of the EWMA chart with fixed control limits,
\code{xtewma.arl} determines the Average Run Length (ARL) by numerically
solving the related ARL integral equation by means of the Nystroem method
based on Gauss-Legendre quadrature.
If \code{limits} is \code{"vacl"}, then the method presented in Knoth (2003) is utilized.
Other values (normal case) for \code{limits} are not yet supported.
}
\value{Except for the fixed limits EWMA charts it returns a single value which resembles the ARL.
For fixed limits charts, it returns a vector of length \code{q} which resembles the ARL and the
sequence of conditional expected delays for \code{q}=1 and \code{q}>1, respectively.}
\references{
K.-H. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving 
average charts, \emph{Appl. Statist. 35}, 151-158.

S. V. Crowder (1987),
A simple method for studying run-length distributions of
exponentially weighted moving average charts,
\emph{Technometrics 29}, 401-407.

J. M. Lucas and M. S. Saccucci (1990),
Exponentially weighted moving average control schemes: Properties
and enhancements, \emph{Technometrics 32}, 1-12.

C. M. Borror, D. C. Montgomery, and G. C. Runger (1999),
Robustness of the EWMA control chart to non-normality ,
\emph{Journal of Quality Technology 31}, 309-316.

S. Knoth (2003),
EWMA schemes with non-homogeneous transition kernels,
\emph{Sequential Analysis 22}, 241-255.

S. Knoth (2004),
Fast initial response features for EWMA Control Charts,
\emph{Statistical Papers 46}, 47-64.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} for zero-state ARL computation of EWMA control charts in the normal case.
}
\examples{
##  Borror/Montgomery/Runger (1999), Table 3
lambda <- 0.1
cE <- 2.703
df <- c(4, 6, 8, 10, 15, 20, 30, 40, 50)
L0 <- rep(NA, length(df))
for ( i in 1:length(df) ) L0[i] <- round(xtewma.arl(lambda, cE*sqrt(df/(df-2)), df, 0), digits=0)
data.frame(df, L0)
}
\keyword{ts}
