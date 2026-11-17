\name{cfa_long_fit_md}
\alias{cfa_long_fit_md}
\title{Fit a longitudinal ordinal confirmatory factor model}
\description{Maximize the observed data likelihood of a longitudinal ordinal CFA model with a probit link.}
\usage{
cfa_long_fit_md(data, namvar, occ = "occ", ID = "ID", Q = NULL, simple = NULL, loadmask = NULL,
    marker = NULL, covariates = NULL, post_covariates = NULL, methInteg = c("QMC",
        "GH", "AGH"), update_gh = 0, nquadqmc = 2000, nquadgh = 8, re_item = TRUE,
    two_tier = FALSE, two_tier_type = c("node_specific", "full_equiv"), theta_time = FALSE,
    meas_invar = NULL, intercepts_free_items = NULL, monotone_tau = TRUE,
    phi_struct = "FREE", maxiter = 500, start = NULL,
    verbose = FALSE, compute_hessian = FALSE, analytic_score = TRUE, hess_eps = 1e-04,
    profile = TRUE, profile_print = TRUE, nquadqmc_re = NULL, re_tier = c("GH",
        "QMC"), parscale = NULL, tol_loglik = 1e-07, tol_par = 0.001, threshold_id = c("tau1_tau2",
        "tau2_tau3"), mean_ident = c("latent_baseline", "marker_intercept"))
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
\item{marker}{One marker item per factor, given as a name or index. Its loading is fixed to one.}
\item{covariates, post_covariates}{Names of numeric predictors of the latent mean and their subset that enters the model after the first occasion.}
\item{methInteg}{Integration method, one of \code{QMC}, \code{GH}, or \code{AGH}. QMC requires \pkg{spacefillr}.}
\item{nquadqmc, nquadqmc_re}{Numbers of Sobol points for the main QMC rule and the optional one dimensional QMC rule for item intercepts.}
\item{re_item}{Include a Gaussian subject by item random intercept.}
\item{two_tier, two_tier_type}{Use an outer rule on the latent trajectory and inner item rules. For AGH, use node-specific adaptation or the full-equivalent block rule.}
\item{re_tier}{Use GH or QMC for the inner rule when the outer rule uses QMC.}
\item{theta_time}{Allow residual variances to vary over time.}
\item{meas_invar}{Four indicators of invariance, in the order thresholds, intercepts, loadings, and residual variances. Time-specific thresholds, intercepts and loadings require the full-equivalent AGH rule.}
\item{intercepts_free_items}{Items whose intercepts are freed after baseline in the full-equivalent model.}
\item{phi_struct}{Use \code{"FREE"} for a free covariance matrix of the stacked latent factors.}
\item{analytic_score, hess_eps}{Use the fixed-grid analytic score where supported, and the finite difference step for Hessian calculation.}
\item{profile, profile_print}{Collect runtime measurements and print their summary.}
\item{threshold_id}{Fix thresholds one and two, or two and three, to zero and one.}
\item{mean_ident}{Fix the baseline factor means to zero or, for full-equivalent AGH, fix the marker intercepts to zero.}
}
\value{An object of class \code{longProbitLVM_fit} with log likelihood, convergence code, parameter estimates and decoded model parameters.}
\details{
The latent covariance is free and the item residuals are conditionally independent. The threshold identification requires at least four categories for each item. The FIML likelihood uses all observed responses within each subject.

For \eqn{T} occasions, \eqn{Q} factors and \eqn{K} included item intercepts, standard GH and AGH integrate \eqn{TQ + K} dimensions jointly. Two-tier integration uses \eqn{TQ} outer dimensions and a one-dimensional rule for each item intercept. Node-specific AGH adapts the inner item rules at each outer node. Full-equivalent AGH uses the joint mode and Hessian, whose dimension is \eqn{TQ + K}. See \code{\link{integration_info}} for the counts used by the fit.
For an invariance sequence based on Wu and Estabrook, see \code{\link{cfa_long_fit_wu2016}}.
}
