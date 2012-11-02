\name{xsewma.arl}
\alias{xsewma.arl}
\title{Compute ARLs of simultaneous EWMA control charts (mean and variance charts)}
\description{Computation of the (zero-state) Average Run Length (ARL)
for different types of simultaneous EWMA control charts 
(based on the sample mean and the sample variance \eqn{S^2})
monitoring normal mean and variance.}
\usage{xsewma.arl(lx, cx, ls, csu, df, mu, sigma, hsx=0, Nx=40, csl=0, hss=1, Ns=40, s2.on=TRUE, sided="upper", qm=30)}
\arguments{
\item{lx}{smoothing parameter lambda of the two-sided mean EWMA chart.}
\item{cx}{control limit of the two-sided mean EWMA control chart.}
\item{ls}{smoothing parameter lambda of the variance EWMA chart.}
\item{csu}{upper control limit of the variance EWMA control chart.}
\item{df}{actual degrees of freedom, corresponds to batch size (for known mean it is equal to the batch size,
for unknown mean it is equal to batch size minus one.}
\item{mu}{true mean.}
\item{sigma}{true standard deviation.}
\item{hsx}{so-called headstart (enables fast initial response) of the mean chart -- 
do not confuse with the true FIR feature considered in xewma.arl; will be updated.}
\item{Nx}{dimension of the approximating matrix of the mean chart.}
\item{csl}{lower control limit of the variance EWMA control chart; default value is 0; not considered if
\code{sided} is \code{"upper"}.}
\item{hss}{headstart (enables fast initial response) of the variance chart.}
\item{Ns}{dimension of the approximating matrix of the variance chart.}
\item{s2.on}{distinguishes between \eqn{S^2}{S^2} and \eqn{S}{S} chart.}
\item{sided}{distinguishes between one- and two-sided two-sided EWMA-\eqn{S^2}{S^2} control charts 
by choosing \code{"upper"} (upper chart without reflection at \code{cl} -- the actual value of
\code{cl} is not used), \code{"Rupper"} (upper chart with reflection at \code{cl}),
\code{"Rlower"} (lower chart with reflection at \code{cu}), and \code{"two"} (two-sided chart), 
respectively.} 
\item{qm}{number of quadrature nodes used for the collocation integrals.}
}
\details{
\code{xsewma.arl} determines the Average Run Length (ARL) by 
an extension of Gan's (derived from ideas already published by Waldmann)
algorithm. The variance EWMA part is treated similarly to the ARL calculation method
deployed for the single variance EWMA charts in Knoth (2005), that is, by means of 
collocation (Chebyshev polynomials). For more details see Knoth (2007).}
\value{Returns a single value which resembles the ARL.}
\references{
K. Waldmann (1986),
Bounds for the distribution of the run length of geometric moving average charts,
\emph{J. R. Stat. Soc., Ser. C, Appl. Stat. 35}, 151-158.

F. F. Gan (1995),
Joint monitoring of process mean and variance using exponentially weighted moving average control charts,
\emph{Technometrics 37}, 446-453.

S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.

S. Knoth (2007),
Accurate ARL calculation for EWMA control charts monitoring simultaneously normal mean and variance,
\emph{Sequential Analysis 26}, 251-264.
}
\author{Sven Knoth}
\seealso{
\code{xewma.arl} and \code{sewma.arl} for zero-state ARL computation of single mean and variance EWMA 
control charts, respectively.
}
\examples{
## Knoth (2007)
## collocation results in Table 1
## Monte Carlo with 10^9 replicates: 252.307 +/- 0.0078

# process parameters
mu <- 0
sigma <- 1
# batch size n=5, df=n-1
df  <- 4
# lambda of mean chart
lx  <- .134
# c_mu^* = .345476571 = cx/sqrt(n) * sqrt(lx/(2-lx)
cx  <- .345476571*sqrt(df+1)/sqrt(lx/(2-lx))
# lambda of variance chart
ls  <- .1
# c_sigma = .477977
csu <- 1 + .477977
# matrix dimensions for mean and variance part
Nx  <- 25
Ns  <- 25
# mode of variance chart
SIDED <- "upper"

arl <- xsewma.arl(lx, cx, ls, csu, df, mu, sigma, Nx=Nx, Ns=Ns, sided=SIDED)
arl
}
\keyword{ts}
