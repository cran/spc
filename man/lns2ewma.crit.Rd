\name{lns2ewma.crit}
\alias{lns2ewma.crit}
\title{Compute critical values of EWMA ln \eqn{S^2}{S^2} control charts (variance charts)}
\description{Computation of the critical values (similar to alarm limits)
for different types of EWMA control charts 
(based on the log of the sample variance \eqn{S^2}) monitoring normal variance.}
\usage{lns2ewma.crit(l,L0,df,sigma0=1,cl=NULL,cu=NULL,hs=NULL,sided="upper",mode="fixed",r=40)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{L0}{in-control ARL.}
\item{df}{actual degrees of freedom, corresponds to subsample size
(for known mean it is equal to the subsample size,
for unknown mean it is equal to subsample size minus one.}
\item{sigma0}{in-control standard deviation.}
\item{cl}{deployed for \code{sided}=\code{"upper"}, that is, upper variance control chart with
lower reflecting barrier \code{cl}.}
\item{cu}{for two-sided (\code{sided}=\code{"two"}) and fixed upper control limit (\code{mode}=\code{"fixed"}),
for all other cases \code{cu} is ignored.}
\item{hs}{so-called headstart (enables fast initial response) -- the default value (hs=NULL) corresponds to the
in-control mean of ln \eqn{S^2}{S^2}.}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA-\eqn{S^2}{S^2} control charts by choosing
\code{"upper"} (upper chart with reflection at \code{cl}), 
\code{"lower"} (lower chart with reflection at \code{cu}), and \code{"two"} (two-sided chart),  respectively.} 
\item{mode}{only deployed for \code{sided}=\code{"two"} -- with \code{"fixed"} an upper control limit (see \code{cu})
is set and only the lower is
calculated to obtain the in-control ARL \code{L0}, while with \code{"unbiased"} a certain unbiasedness of the ARL
function is guaranteed (here, both the
lower and the upper control limit are calculated). With \code{"vanilla"} limits symmetric around the in-control
mean of ln \eqn{S^2}{S^2}
are determined, while for \code{"eq.tails"} the in-control ARL values of two single EWMA variance charts
(decompose the two-sided scheme into one lower and one upper scheme) are matched.}
\item{r}{dimension of the resulting linear equation system: the larger the more accurate.}
}
\details{
\code{lns2ewma.crit} determines the critical values (similar to alarm limits) for given in-control ARL \code{L0}
by applying secant rule and using \code{lns2ewma.arl()}.
In case of \code{sided}=\code{"two"} and \code{mode}=\code{"unbiased"}
a two-dimensional secant rule is applied that also ensures that the
maximum of the ARL function for given standard deviation is attained
at \code{sigma0}. See Knoth (2010) and the related example.
}
\value{Returns the lower and upper control limit \code{cl} and \code{cu}.}
\references{
C. A. Acosta-Mej\'ia and J. J. Pignatiello Jr. and B. V. Rao (1999),
A comparison of control charting procedures for monitoring process dispersion,
\emph{IIE Transactions 31}, 569-579.

S. V. Crowder and M. D. Hamilton (1992),
An EWMA for monitoring a process standard deviation,
\emph{Journal of Quality Technology 24}, 12-21.

S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.

S. Knoth (2010),
Control Charting Normal Variance -- Reflections, Curiosities, and Recommendations,
in \emph{Frontiers in Statistical Quality Control 9},
H.-J. Lenz and P.-T. Wilrich (Eds.),
Physica Verlag, Heidelberg, Germany, 3-18.
}
\author{Sven Knoth}
\seealso{\code{lns2ewma.arl} for calculation of ARL of EWMA ln \eqn{S^2}{S^2} control charts.}
\examples{
## Knoth (2005)
## compare with 1.05521 mentioned on page 350 third line from below
L0 <- 200
lambda <- .05
df <- 4
limits <- lns2ewma.crit(lambda, L0, df, cl=0, hs=0)
limits["cu"]/sqrt( lambda/(2-lambda)*(2/df+2/df^2+4/3/df^3-16/15/df^5) )
}
\keyword{ts}
