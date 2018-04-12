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
could be found at
%\url{http://lib.stat.cmu.edu/jqt/33-4}.
http://lib.stat.cmu.edu/jqt/33-4.
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
\emph{Journal of Quality Technology 31(1)}, 120-123.

Giovanna Capizzi and Guido Masarotto (2010),
Evaluation of the run-length distribution for a combined Shewhart-EWMA control chart,
\emph{Statistics and Computing 20(1)}, 23-33.

Sven Knoth (2017),
ARL Numerics for MEWMA Charts,
\emph{Journal of Quality Technology 49(1)}, 78-89.

Wade E. Molnau et al. (2001),
A Program for ARL Calculation for Multivariate EWMA Charts,
\emph{Journal of Quality Technology 33(4)}, 515-521.

Sharad S. Prabhu and George C. Runger (1997),
Designing a multivariate EWMA control chart,
\emph{Journal of Quality Technology 29(1)}, 8-15.

Steven E. Rigdon (1995a), An integral equation for the in-control average run length of a multivariate
exponentially weighted moving average control chart, \emph{J. Stat. Comput. Simulation 52(4)}, 351-365.

Steven E. Rigdon (1995b), A double-integral equation for the average run length of a multivariate
exponentially weighted moving average control chart, \emph{Stat. Probab. Lett. 24(4)}, 365-373.

George C. Runger and Sharad S. Prabhu (1996),
A Markov Chain Model for the Multivariate Exponentially Weighted Moving Averages Control Chart,
\emph{J. Amer. Statist. Assoc. 91(436)}, 1701-1706.
}
\author{Sven Knoth}
\seealso{
\code{mewma.crit} for getting the alarm threshold to attain a certain in-control ARL.
}
\examples{
# Rigdon (1995a), p. 357, Tab. 1
p <- 2
r <- 0.25
h4 <- c(8.37, 9.90, 11.89, 13.36, 14.82, 16.72)
for ( i in 1:length(h4) ) cat(paste(h4[i], "\t", round(mewma.arl(r, h4[i], p, ntype="ra")), "\n"))

r <- 0.1
h4 <- c(6.98, 8.63, 10.77, 12.37, 13.88, 15.88)
for ( i in 1:length(h4) ) cat(paste(h4[i], "\t", round(mewma.arl(r, h4[i], p, ntype="ra")), "\n"))


# Rigdon (1995b), p. 372, Tab. 1
\dontrun{
r <- 0.1
p <- 4
h <- 12.73
for ( sdelta in c(0, 0.125, 0.25, .5, 1, 2, 3) )
  cat(paste(sdelta, "\t",
      round(mewma.arl(r, h, p, delta=sdelta^2, ntype="ra", r=25), digits=2), "\n"))

p <- 5
h <- 14.56
for ( sdelta in c(0, 0.125, 0.25, .5, 1, 2, 3) )
  cat(paste(sdelta, "\t",
      round(mewma.arl(r, h, p, delta=sdelta^2, ntype="ra", r=25), digits=2), "\n"))

p <- 10
h <- 22.67
for ( sdelta in c(0, 0.125, 0.25, .5, 1, 2, 3) )
  cat(paste(sdelta, "\t",
      round(mewma.arl(r, h, p, delta=sdelta^2, ntype="ra", r=25), digits=2), "\n"))
}

# Runger/Prabhu (1996), p. 1704, Tab. 1
\dontrun{
r <- 0.1
p <- 4
H <- 12.73
cat(paste(0, "\t", round(mewma.arl(r, H, p, delta=0, ntype="mc", r=50), digits=2), "\n"))
for ( delta in c(.5, 1, 1.5, 2, 3) )
  cat(paste(delta, "\t",
      round(mewma.arl(r, H, p, delta=delta, ntype="mc", r=25), digits=2), "\n"))
# compare with Fortran program (MEWMA-ARLs.f90) from Molnau et al. (2001) with m1 = m2 = 25
# H4      P     R   DEL  ARL
# 12.73  4.  0.10  0.00 199.78
# 12.73  4.  0.10  0.50  35.05
# 12.73  4.  0.10  1.00  12.17
# 12.73  4.  0.10  1.50   7.22
# 12.73  4.  0.10  2.00   5.19
# 12.73  4.  0.10  3.00   3.42

p <- 20
H <- 37.01
cat(paste(0, "\t",
    round(mewma.arl(r, H, p, delta=0, ntype="mc", r=50), digits=2), "\n"))
for ( delta in c(.5, 1, 1.5, 2, 3) )
  cat(paste(delta, "\t",
      round(mewma.arl(r, H, p, delta=delta, ntype="mc", r=25), digits=2), "\n"))
# compare with Fortran program (MEWMA-ARLs.f90) from Molnau et al. (2001) with m1 = m2 = 25
# H4      P     R   DEL  ARL
# 37.01 20.  0.10  0.00 199.09
# 37.01 20.  0.10  0.50  61.62
# 37.01 20.  0.10  1.00  20.17
# 37.01 20.  0.10  1.50  11.40
# 37.01 20.  0.10  2.00   8.03
# 37.01 20.  0.10  3.00   5.18
}

# Knoth (2017), p. 85, Tab. 3, rows with p=3
\dontrun{
p <- 3
lambda <- 0.05
h4 <- mewma.crit(lambda, 200, p)
benchmark <- mewma.arl(lambda, h4, p, delta=1, r=50)
  
mc.arl  <- mewma.arl(lambda, h4, p, delta=1, r=25, ntype="mc")
ra.arl  <- mewma.arl(lambda, h4, p, delta=1, r=27, ntype="ra")
co.arl  <- mewma.arl(lambda, h4, p, delta=1, r=12, ntype="co2")
gl3.arl <- mewma.arl(lambda, h4, p, delta=1, r=30, ntype="gl3")
gl5.arl <- mewma.arl(lambda, h4, p, delta=1, r=25, ntype="gl5")
  
abs( benchmark - data.frame(mc.arl, ra.arl, co.arl, gl3.arl, gl5.arl) )
}

# Prabhu/Runger (1997), p. 13, Tab. 3
\dontrun{
p <- 2
r <- 0.1
H <- 8.64
cat(paste(0, "\t",
    round(mewma.ad(r, H, p, delta=0, type="cycl", ntype="mc", r=60), digits=2), "\n"))
for ( delta in c(.5, 1, 1.5, 2, 3) )
  cat(paste(delta, "\t",
      round(mewma.ad(r, H, p, delta=delta, type="cycl", ntype="mc", r=30), digits=2), "\n"))

# better accuracy
for ( delta in c(0, .5, 1, 1.5, 2, 3) )
  cat(paste(delta, "\t",
      round(mewma.ad(r, H, p, delta=delta^2, type="cycl", r=30), digits=2), "\n"))
}
}
\keyword{ts}
