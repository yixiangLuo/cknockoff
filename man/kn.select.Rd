% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knockoff.R
\name{kn.select}
\alias{kn.select}
\title{Select variables by knockoff.}
\usage{
kn.select(kn_statistics, alpha, selective = T, early_stop = F)
}
\arguments{
\item{kn_statistics}{The knockoff feature statistics (W-statistics)}

\item{alpha}{target false discovery rate}

\item{selective}{If TRUE, use selective SeqStep for selection
(only select variables with a positive W-statistics);
if FALSE, use SeqStep for selection
(select variables if their W-statistics have large absolute values.)}

\item{early_stop}{If FALSE, the selection set is the same as in knockoff
(stop at \eqn{\tau});
If TRUE, the selection set is what we used to construct the budget \eqn{b_j}
(stop at \eqn{\tau_1}).}
}
\value{
A list containing the following components:
\item{selected}{vector of selected variables.}
\item{fdp_est}{the estimated FDP.}
\item{W_k_hat}{The absolute value of the W-statistics where we stop.}
}
\description{
Select variables relevant for predicting the outcome of interest, using the
knockoff feature statistics and applying selective SeqStep.
}
\examples{

}
\keyword{internal}
