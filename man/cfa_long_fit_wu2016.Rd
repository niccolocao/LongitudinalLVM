\name{cfa_long_fit_wu2016}
\alias{cfa_long_fit_wu2016}
\alias{measurement_invariance_wu2016}
\title{Fit longitudinal ordinal measurement invariance models}
\description{Fit the ordinal CFA invariance restrictions based on Wu and Estabrook (2016) using the full-equivalent adaptive quadrature rule.}
\usage{
cfa_long_fit_wu2016(data, namvar, occ = "occ", ID = "ID", Q = NULL, simple = NULL, loadmask = NULL,
    stage = c("baseline", "threshold", "threshold_loading", "threshold_loading_intercept",
        "strict"), nquadgh = 7L, update_gh = 0L, maxiter = 500L, monotone_tau = TRUE,
    start = NULL, parscale = NULL, tol_loglik = 1e-07, tol_par = 0.001, compute_hessian = FALSE,
    verbose = FALSE)
measurement_invariance_wu2016(data, namvar, occ = "occ", ID = "ID", Q = NULL, simple = NULL, loadmask = NULL,
    stages = c("baseline", "threshold", "threshold_loading", "threshold_loading_intercept",
        "strict"), ...)
}
\arguments{
\item{data}{A data frame with one row per subject and occasion. Item responses can be integer categories starting at 0 or 1, or ordered factors. Missing responses are coded as \code{NA}.}
\item{namvar}{Names of the item columns.}
\item{occ, ID}{Names of the occasion and subject columns. Occasions are ordered by their labels.}
\item{Q, simple, loadmask}{Number of factors, one factor index per item, or a matrix of active loadings with one row per item and one column per factor. The default is one factor with all items loading on it.}
\item{nquadgh}{Number of Gauss-Hermite nodes per dimension.}
\item{update_gh}{Use zero to update the adaptive grid at each objective evaluation. A positive value sets the initial number of iterations per fixed grid block.}
\item{maxiter}{Iteration budget. Zero evaluates the likelihood at the starting values.}
\item{monotone_tau}{Use the logarithms of positive increments to keep thresholds ordered.}
\item{start}{Optional starting vector on the optimization scale.}
\item{parscale}{Optional parameter scaling for BFGS.}
\item{tol_loglik, tol_par}{Relative convergence tolerances for the objective and parameters.}
\item{compute_hessian}{Compute the observed information at the reported fit.}
\item{verbose}{Print optimization progress.}
\item{stage, stages}{One invariance stage or a vector of stages in order from weaker to stronger restrictions.}
\item{...}{Further arguments passed to cfa_long_fit_wu2016.}
}
\value{A single fitted model or an object containing the fitted models and a table of adjacent likelihood ratio comparisons.}
\details{
Each item must load on exactly one factor and have at least three categories.
A comparison p-value is provided for converged fits with a positive difference in degrees of freedom and nonnegative likelihood ratio to numerical tolerance.
}
