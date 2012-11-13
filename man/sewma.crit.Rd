\name{sewma.crit}
\alias{sewma.crit}
\title{Compute critical values of EWMA control charts (variance charts)}
\description{Computation of the critical values (similar to alarm limits)
for different types of EWMA control charts (based on the sample variance
\eqn{S^2}) monitoring normal variance.}
\usage{sewma.crit(l,L0,df,sigma0=1,cl=NULL,cu=NULL,hs=NULL,s2.on=TRUE,
sided="upper",mode="fixed",ur=4,r=40,qm=30)}
\arguments{
\item{l}{smoothing parameter lambda of the EWMA control chart.}
\item{L0}{in-control ARL.}
\item{df}{actual degrees of freedom, corresponds to subgroup size (for known mean it is equal to the subgroup size,
for unknown mean it is equal to subgroup size minus one.}
\item{sigma0}{in-control standard deviation.}
\item{cl}{deployed for \code{sided}=\code{"Rupper"}, that is, upper variance control chart
with lower reflecting barrier \code{cl}.}
\item{cu}{for two-sided (\code{sided}=\code{"two"}) and fixed upper
control limit (\code{mode}=\code{"fixed"}) a value larger than \code{sigma0}
has to been given, for all other cases \code{cu} is ignored.}
\item{hs}{so-called headstart (enables fast initial response); the default (\code{NULL})
yields the expected in-control value of \eqn{S^2}{S^2} (1) and \eqn{S}{S} (\eqn{c_4}{c_4}),
respectively.}
\item{s2.on}{distinguishes between \eqn{S^2}{S^2} and \eqn{S}{S} chart.}
\item{sided}{distinguishes between one- and two-sided
two-sided EWMA-\eqn{S^2}{S^2} control charts
by choosing \code{"upper"} (upper chart
without reflection at \code{cl} -- the actual value of
\code{cl} is not used),
\code{"Rupper"} (upper chart with reflection at \code{cl}),
\code{"Rlower"} (lower chart with reflection at \code{cu}),
and \code{"two"} (two-sided chart),
respectively.}
\item{mode}{only deployed for \code{sided}=\code{"two"} -- with \code{"fixed"}
an upper control limit (see \code{cu}) is set and only the lower is
calculated to obtain the in-control ARL \code{L0}, while with \code{"unbiased"}
a certain unbiasedness of the ARL function is guaranteed (here, both the
lower and the upper control limit are calculated). With \code{"vanilla"} limits symmetric
around 1 (the in-control value of the variance)
are determined, while for \code{"eq.tails"} the in-control ARL values of two single EWMA
variance charts (decompose the two-sided
scheme into one lower and one upper scheme) are matched.}
\item{ur}{truncation of lower chart for \code{eq.tails} mode.}
\item{r}{dimension of the resulting linear equation system
(highest order of the collocation polynomials).}
\item{qm}{number of quadrature nodes for calculating the collocation definite integrals.}
}
\details{
\code{sewma.crit} determines the critical values (similar to alarm limits)
for given in-control ARL \code{L0}
by applying secant rule and using \code{sewma.arl()}.
In case of \code{sided}=\code{"two"} and \code{mode}=\code{"unbiased"}
a two-dimensional secant rule is applied that also ensures that the
maximum of the ARL function for given standard deviation is attained
at \code{sigma0}. See Knoth (2010) and the related example.
}
\value{Returns the lower and upper control limit \code{cl} and \code{cu}.}
\references{
H.-J. Mittag and D. Stemann and B. Tewes (1998),
EWMA-Karten zur \"Uberwachung der Streuung von Qualit\"atsmerkmalen,
\emph{Allgemeines Statistisches Archiv 82}, 327-338,

C. A. Acosta-Mej\'ia and J. J. Pignatiello Jr. and B. V. Rao (1999),
A comparison of control charting procedures for monitoring process dispersion,
\emph{IIE Transactions 31}, 569-579.

S. Knoth (2005),
Accurate ARL computation for EWMA-\eqn{S^2}{S^2} control charts,
\emph{Statistics and Computing 15}, 341-352.

S. Knoth (2006a),
Computation of the ARL for CUSUM-\eqn{S^2}{S^2} schemes,
\emph{Computational Statistics & Data Analysis 51}, 499-512.

S. Knoth (2006b), 
The art of evaluating monitoring schemes -- how to measure the performance of control charts?
in \emph{Frontiers in Statistical Quality Control 8},
H.-J. Lenz and P.-T. Wilrich (Eds.),
Physica Verlag, Heidelberg, Germany, 74-99.

S. Knoth (2010),
Control Charting Normal Variance -- Reflections, Curiosities, and Recommendations,
in \emph{Frontiers in Statistical Quality Control 9},
H.-J. Lenz and P.-T. Wilrich (Eds.),
Physica Verlag, Heidelberg, Germany, 3-18.
}
\author{Sven Knoth}
\seealso{\code{sewma.arl} for calculation of ARL of variance charts.}
\examples{
## Mittag et al. (1998)
## compare their upper critical value 2.91 that
## leads to the upper control limit via the formula shown below
## (for the usual upper EWMA \eqn{S^2}{S^2}).
## See Knoth (2006b) for a discussion of this EWMA setup and it's evaluation.

l  <- 0.18
L0 <- 250
df <- 4
limits <- sewma.crit(l, L0, df)
limits["cu"]

limits.cu.mittag_et_al <- 1 + sqrt(l/(2-l))*sqrt(2/df)*2.91
limits.cu.mittag_et_al

## Knoth (2005)
## reproduce the critical value given in Figure 2 (c=1.661865) for
## upper EWMA \eqn{S^2}{S^2} with df=1

l  <- 0.025
L0 <- 250
df <- 1
limits <- sewma.crit(l, L0, df)
cv.Fig2 <- (limits["cu"]-1)/( sqrt(l/(2-l))*sqrt(2/df) )
cv.Fig2

## the small difference (sixth digit after decimal point) stems from
## tighter criterion in the secant rule implemented in the R package.

## demo of unbiased ARL curves
## Deploy, please, not matrix dimensions smaller than 50 -- for the
## sake of accuracy, the value 80 was used.
## Additionally, this example needs between 1 and 2 minutes on a 1.6 Ghz box. 

\dontrun{
l  <- 0.1
L0 <- 500
df <- 4
limits <- sewma.crit(l, L0, df, sided="two", mode="unbiased", r=80)
SEWMA.arl <- Vectorize(sewma.arl, "sigma")
SEWMA.ARL <- function(sigma) 
  SEWMA.arl(l, limits[1], limits[2], sigma, df, sided="two", r=80)
layout(matrix(1:2, nrow=1))
curve(SEWMA.ARL, .75, 1.25, log="y")
curve(SEWMA.ARL, .95, 1.05, log="y")}
# the above stuff needs about 1 minute

## control limits for upper and lower EWMA charts with reflecting barriers
## (reflection at in-control level sigma0 = 1)
## examples from Knoth (2006a), Tables 4 and 5

\dontrun{
## upper chart with reflection at sigma0=1 in Table 4: c = 2.4831
l <- 0.15
L0 <- 100
df <- 4
limits <- sewma.crit(l, L0, df, cl=1, sided="Rupper", r=100)
cv.Tab4 <- (limits["cu"]-1)/( sqrt(l/(2-l))*sqrt(2/df) )
cv.Tab4

## lower chart with reflection at sigma0=1 in Table 5: c = 2.0613
l <- 0.115
L0 <- 200
df <- 5
limits <- sewma.crit(l, L0, df, cu=1, sided="Rlower", r=100)
cv.Tab5 <- -(limits["cl"]-1)/( sqrt(l/(2-l))*sqrt(2/df) )
cv.Tab5}
}
\keyword{ts}
