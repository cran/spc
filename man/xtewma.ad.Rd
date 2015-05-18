\name{xtewma.ad}
\alias{xtewma.ad}
\title{Compute steady-state ARLs of EWMA control charts, t distributed data}
\description{Computation of the steady-state Average Run Length (ARL)
for different types of EWMA control charts monitoring the mean of t distributed data.}
\usage{xtewma.ad(l, c, df, mu1, mu0=0, zr=0, z0=0, sided="one", limits="fix",
steady.state.mode="conditional", mode="tan", r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{c}{critical value (similar to alarm limit) of the EWMA control chart.}
\item{df}{degrees of freedom -- parameter of the t distribution.}
\item{mu1}{in-control mean.}
\item{mu0}{out-of-control mean.}
\item{zr}{reflection border for the one-sided chart.}
\item{z0}{restarting value of the EWMA sequence in case of a false alarm in
\code{steady.state.mode="cyclical"}.}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA control
chart by choosing \code{"one"} and \code{"two"}, respectively.}
\item{limits}{distinguishes between different control limits behavior.}
\item{steady.state.mode}{distinguishes between two steady-state modes -- conditional and cyclical.}
\item{mode}{Controls the type of variables substitution that might improve the numerical performance. Currently,
\code{"identity"}, \code{"sin"}, \code{"sinh"}, and \code{"tan"} (default) are provided.}
\item{r}{number of quadrature nodes, dimension of the resulting linear
equation system is equal to \code{r+1} (one-sided) or \code{r}
(two-sided).}
}
\details{
\code{xewma.ad} determines the steady-state Average Run Length (ARL)
by numerically solving the related ARL integral equation by means
of the Nystroem method based on Gauss-Legendre quadrature
and using the power method for deriving the largest in magnitude
eigenvalue and the related left eigenfunction.
}
\value{Returns a single value which resembles the steady-state ARL.}
\references{
R. B. Crosier (1986),
A new two-sided cumulative quality control scheme,
\emph{Technometrics 28}, 187-194.

S. V. Crowder (1987),
A simple method for studying run-length distributions of exponentially weighted
moving average charts,
\emph{Technometrics 29}, 401-407.

J. M. Lucas and M. S. Saccucci (1990),
Exponentially weighted moving average control schemes: Properties and enhancements,
\emph{Technometrics 32}, 1-12.
}
\author{Sven Knoth}
\seealso{
\code{xtewma.arl} for zero-state ARL computation and
\code{xewma.ad} for the steady-state ARL for normal data.}
\examples{
## will follow
}
\keyword{ts}
