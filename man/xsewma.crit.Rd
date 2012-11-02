\name{xsewma.crit}
\alias{xsewma.crit}
\title{Compute critical values of simultaneous EWMA control charts (mean and variance charts)}
\description{Computation of the critical values (similar to alarm limits)
for different types of simultaneous EWMA control charts 
(based on the sample mean and the sample variance \eqn{S^2})
monitoring normal mean and variance.}
\usage{xsewma.crit(lx, ls, L0, df, mu0=0, sigma0=1, cu=NULL, hsx=0, hss=1, s2.on=TRUE, sided="upper", mode="fixed", Nx=20, Ns=40, qm=30)}
\arguments{
\item{lx}{smoothing parameter lambda of the two-sided mean EWMA chart.}
\item{ls}{smoothing parameter lambda of the variance EWMA chart.}
\item{L0}{in-control ARL.}
\item{mu0}{in-control mean.}
\item{sigma0}{in-control standard deviation.}
\item{cu}{for two-sided (\code{sided}=\code{"two"}) and fixed upper
control limit (\code{mode}=\code{"fixed"}) a value larger than \code{sigma0}
has to been given, for all other cases \code{cu} is ignored.}
\item{hsx}{so-called headstart (enables fast initial response) of the mean chart -- 
do not confuse with the true FIR feature considered in xewma.arl; will be updated.}
\item{hss}{headstart (enables fast initial response) of the variance chart.}
\item{df}{actual degrees of freedom, corresponds to batch size (for known mean it is equal to the batch size,
for unknown mean it is equal to batch size minus one.}
\item{s2.on}{distinguishes between \eqn{S^2}{S^2} and \eqn{S}{S} chart.}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA-\eqn{S^2}{S^2} control charts
by choosing \code{"upper"} (upper chart without reflection at \code{cl} -- the actual value of
\code{cl} is not used), \code{"Rupper"} (upper chart with reflection at \code{cl}),
\code{"Rlower"} (lower chart with reflection at \code{cu}), and \code{"two"} (two-sided chart),
respectively.}
\item{mode}{only deployed for \code{sided}=\code{"two"} -- with \code{"fixed"}
an upper control limit (see \code{cu}) is set and only the lower is
determined to obtain the in-control ARL \code{L0}, while with \code{"unbiased"}
a certain unbiasedness of the ARL function is guaranteed (here, both the
lower and the upper control limit are calculated).}
\item{Nx}{dimension of the approximating matrix of the mean chart.}
\item{Ns}{dimension of the approximating matrix of the variance chart.}
\item{qm}{number of quadrature nodes used for the collocation integrals.}
}
\details{
\code{xsewma.crit} determines the critical values (similar to alarm limits)
for given in-control ARL \code{L0}
by applying secant rule and using \code{xsewma.arl()}.
In case of \code{sided}=\code{"two"} and \code{mode}=\code{"unbiased"}
a two-dimensional secant rule is applied that also ensures that the
maximum of the ARL function for given standard deviation is attained
at \code{sigma0}. See Knoth (2007) for details and application.
}
\value{Returns the critical value of the two-sided mean EWMA chart and
the lower and upper controls limit \code{cl} and \code{cu} of the variance EWMA chart.}
\references{

S. Knoth (2007),
Accurate ARL calculation for EWMA control charts monitoring simultaneously normal mean and variance,
\emph{Sequential Analysis 26}, 251-264.
}
\author{Sven Knoth}
\seealso{\code{xsewma.arl} for calculation of ARL of simultaneous EWMA charts.}
\examples{
## Knoth (2007)
## results in Table 2

# batch size n=5, df=n-1
df  <- 4
# lambda of mean chart
lx  <- .134
# lambda of variance chart
ls  <- .1
# in-control ARL 
L0 <- 252.3
# matrix dimensions for mean and variance part
Nx  <- 25
Ns  <- 25
# mode of variance chart
SIDED <- "upper"

crit <- xsewma.crit(lx, ls, L0, df, sided=SIDED, Nx=Nx, Ns=Ns)
crit

## output as used in Knoth (2007)
crit["cx"]/sqrt(df+1)*sqrt(lx/(2-lx))
crit["cu"] - 1
}
\keyword{ts}
