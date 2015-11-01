\name{mewma.arl}
\alias{mewma.arl}
\alias{mewma.arl.f}
\alias{mewma.ad}
\title{Compute ARLs of MEWMA control charts}
\description{Computation of the (zero-state) Average Run Length (ARL)
for multivariate exponentially weighted moving average (MEWMA) charts monitoring multivariate normal mean.}
\usage{mewma.arl(l, cE, p, delta=0, hs=0, r=20, ntype=NULL, qm0=20, qm1=qm0)

mewma.arl.f(l, cE, p, delta=0, r=20, ntype=NULL, qm0=20, qm1=qm0)

mewma.ad(l, cE, p, delta=0, r=20, n=20, type="cond", hs=0, ntype=NULL, qm0=20, qm1=qm0)}
\arguments{
\item{l}{smoothing parameter lambda of the MEWMA control chart.}
\item{cE}{alarm threshold of the MEWMA control chart.}
\item{p}{dimension of multivariate normal distribution.}
\item{delta}{magnitude of the potential change, \code{delta=0} refers to the in-control state.}
\item{hs}{so-called headstart (enables fast initial response) -- must be non-negative.}
\item{r}{number of quadrature nodes -- dimension of the resulting linear equation system
for \code{delta} = 0. For non-zero \code{delta} this dimension is mostly r^2 (Markov chain approximation leads
to some larger values). Caution: If \code{ntype} is set to \code{"co"} (collocation), then values of \code{r}
larger than 20 lead to large computing times.
For the other selections this would happen for values larger than 40.}
\item{ntype}{choose the numerical algorithm to solve the ARL integral equation. For \code{delta}=0:
Possible values are
\code{"gl"}, \code{"gl2"} (gauss-legendre, classic and with variables change: square),
\code{"co"} (collocation, for \code{delta} > 0 with sin transformation),
\code{"ra"} (radau),
\code{"cc"} (clenshaw-curtis),
\code{"mc"} (markov chain),
and \code{"sr"} (simpson rule).
For \code{delta} larger than 0, some more values besides the others are possible:
\code{"gl3"}, \code{"gl4"}, \code{"gl5"} (gauss-legendre with a further change in variables: sin, tan, sinh),
\code{"co2"}, \code{"co3"} (collocation with some trimming and tan as quadrature stabilizing transformations, respectively).
If it is set to \code{NULL} (the default), then for \code{delta}=0 then \code{"gl2"} is chosen.
If \code{delta} larger than 0, then for \code{p} equal 2 or 4 \code{"gl3"} and for all other values \code{"gl5"} is taken.
\code{"ra"} denotes the method used in Rigdon (1995a). \code{"mc"} denotes the Markov chain approximation.}
\item{type}{switch between \code{"cond"} and \code{"cycl"} for differentiating between the conditional
(no false alarm) and the cyclical (after false alarm re-start in \code{hs}), respectively.}
\item{n}{number of quadrature nodes for Calculating the steady-state ARL integral(s).}
\item{qm0,qm1}{number of collocation quadrature nodes for the out-of-control case (\code{qm0} for the inner integral,
\code{qm1} for the outer one), that is, for positive \code{delta},
and for the in-control case (now only \code{qm0} is deployed) if via \code{ntype} the collocation procedure is requested.}
}
\details{Basically, this is the implementation of different numerical algorithms for
solving the integral equation for the MEWMA in-control (\code{delta} = 0) ARL introduced in Rigdon (1995a)
and out-of-control (\code{delta} != 0) ARL in Rigdon (1995b).
Most of them are nothing else than the Nystroem approach -- the integral is replaced by a suitable quadrature.
Here, the Gauss-Legendre (more powerful), Radau (used by Rigdon, 1995a), Clenshaw-Curtis, and
Simpson rule (which is really bad) are provided.
Additionally, the collocation approach is offered as well, because it is much better for small odd values for \code{p}.
FORTRAN code for the Radau quadrature based Nystroem of Rigdon (1995a)
was published in Bodden and Rigdon (1999) -- see also \url{http://lib.stat.cmu.edu/jqt/31-1}.
Furthermore, FORTRAN code for the Markov chain approximation (in- and out-ot-control)
could be found at \url{http://lib.stat.cmu.edu/jqt/33-4}.
The related papers are Runger and Prabhu (1996) and Molnau et al. (2001).
The idea of the Clenshaw-Curtis quadrature was taken from
Capizzi and Masarotto (2010), who successfully deployed a modified Clenshaw-Curtis quadrature
to calculate the ARL of combined (univariate) Shewhart-EWMA charts. It turns out that it works also nicely for the
MEWMA ARL. The version \code{mewma.arl.f()} without the argument \code{hs} provides the ARL as function of one (in-control)
or two (out-of-control) arguments.
}
\value{Returns a single value which is simply the zero-state ARL.}
\references{
Kevin M. Bodden and Steven E. Rigdon (1999),
A program for approximating the in-control ARL for the MEWMA chart,
\emph{Journal of Quality Technology 31}, 120-123.

Giovanna Capizzi and Guido Masarotto (2010),
Evaluation of the run-length distribution for a combined Shewhart-EWMA control chart,
\emph{Statistics and Computing 20}, 23-33.

Wade E. Molnau et al. (2001),
A Program for ARL Calculation for Multivariate EWMA Charts,
\emph{Journal of Quality Technology 33}, 515-521.

Steven E. Rigdon (1995a), An integral equation for the in-control average run length of a multivariate
exponentially weighted moving average control chart, \emph{J. Stat. Comput. Simulation 52}, 351-365.

Steven E. Rigdon (1995b), A double-integral equation for the average run length of a multivariate
exponentially weighted moving average control chart, \emph{Stat. Probab. Lett. 24}, 365-373.

George C. Runger and Sharad S. Prabhu (1996),
A Markov Chain Model for the Multivariate Exponentially Weighted Moving Averages Control Chart,
\emph{J. Amer. Statist. Assoc. 91}, 1701-1706.
}
\author{Sven Knoth}
\seealso{
\code{mewma.crit} for getting the alarm threshold to attain a certain in-control ARL.
}
\examples{
# Rigdon (1995b), p. 372, Tab. 1
r <- 0.1
p <- 4
h <- 12.73
L0a <- mewma.arl(r, h, p) # defaults to "gl2" because of the even p.
L0b <- mewma.arl(r, h, p, ntype="co")
L0c <- mewma.arl(r, h, p, r=48, ntype="ra")
L0d <- mewma.arl(r, h, p, ntype="cc")
L0e <- mewma.arl(r, h, p, ntype="mc")
data.frame(L0a, L0b, L0c, L0d, L0e)

# original (Rigdon 1995a) implicite value is 200
p <- 3
h <- 14.98
L0a  <- mewma.arl(r, h, p, ntype="gl2")
L0aa <- mewma.arl(r, h, p, r=48, ntype="gl2")
L0b  <- mewma.arl(r, h, p, ntype="co")
L0bb <- mewma.arl(r, h, p, r=48, ntype="co")
L0c  <- mewma.arl(r, h, p, r=48, ntype="ra")
L0d  <- mewma.arl(r, h, p, ntype="cc")
L0dd <- mewma.arl(r, h, p, r=48, ntype="cc")
L0e  <- mewma.arl(r, h, p, ntype="mc")
L0ee <- mewma.arl(r, h, p, r=48, ntype="mc")
data.frame(L0a, L0aa, L0b, L0bb, L0c, L0d, L0dd, L0e, L0ee)
# original (Rigdon 1995a) implicite value is 1000

# Rigdon (1995b), p. 372, Tab. 1
p <- 5
h <- 14.56
L1 <- mewma.arl(r, h, p, delta=1, r=20)
L1
# original value is 12.9
}
\keyword{ts}
