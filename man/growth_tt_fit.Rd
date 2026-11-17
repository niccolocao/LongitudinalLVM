\name{growth_tt_fit}
\alias{growth_tt_fit}
\alias{growth_tt_fit_flexible}
\title{Fit a general longitudinal ordinal latent growth model}
\description{Maximize the observed data likelihood for an ordinal probit latent growth model.}
\usage{
growth_tt_fit(data, namvar, occ = "occ", ID = "ID", loadmask = NULL, factor_names = NULL,
    identification = c("effects", "baseline_std", "marker", "mixed_marker"),
    marker_item = NULL, threshold_identification = c("first_two", "central_01"),
    intercept_time_free = NULL, growth_shape = c("offset", "rho"), rho_fixed = NULL,
    methInteg = c("GH", "QMC", "AGH", "AGH_FULL_EQUIV", "AGH_RE_NESTED",
        "AGH_RE_FULL_EQUIV"), nquadgh = 8, nquadqmc = 2000, update_gh = 1L,
    analytic_score = TRUE, re_item = TRUE, theta_time = TRUE, delta_free = TRUE,
    monotone_tau = FALSE, simple_growth = TRUE, same_design = TRUE, subject_specific_design = NULL,
    time_scores = NULL, fixed_covariates = NULL, time_interactions = NULL,
    intercept_covariates = NULL, slope_covariates = NULL, X = NULL, Z = NULL,
    delta_mask = NULL, fixed_names = NULL, random_names = NULL, maxiter = 500,
    start = NULL, verbose = FALSE, compute_hessian = FALSE, hessian_method = NULL,
    hess_eps = 1e-04, profile = TRUE, profile_print = TRUE, parscale = NULL,
    debug_par = FALSE, debug_par_each_eval = TRUE, debug_par_full = TRUE,
    debug_par_max_eval = Inf, tol_loglik = 1e-07, tol_par = 0.001, Q = NULL,
    simple = NULL, nquadgh_u = nquadgh, nquadgh_zeta = nquadgh, nquadgh_re = nquadgh,
    rho_pattern = NULL, two_tier = NULL)
growth_tt_fit_flexible(data, namvar, occ = "occ", ID = "ID", loadmask = NULL, factor_names = NULL,
    identification = c("effects", "baseline_std", "marker", "mixed_marker"),
    marker_item = NULL, threshold_identification = c("first_two", "central_01"),
    intercept_time_free = NULL, growth_shape = c("offset", "rho"), rho_fixed = NULL,
    methInteg = c("GH", "QMC", "AGH", "AGH_FULL_EQUIV", "AGH_RE_NESTED",
        "AGH_RE_FULL_EQUIV"), nquadgh = 8, nquadqmc = 2000, update_gh = 1L,
    analytic_score = TRUE, re_item = TRUE, theta_time = TRUE, delta_free = TRUE,
    monotone_tau = FALSE, simple_growth = TRUE, same_design = TRUE, subject_specific_design = NULL,
    time_scores = NULL, fixed_covariates = NULL, time_interactions = NULL,
    intercept_covariates = NULL, slope_covariates = NULL, X = NULL, Z = NULL,
    delta_mask = NULL, fixed_names = NULL, random_names = NULL, maxiter = 500,
    start = NULL, verbose = FALSE, compute_hessian = FALSE, hessian_method = NULL,
    hess_eps = 1e-04, profile = TRUE, profile_print = TRUE, parscale = NULL,
    debug_par = FALSE, debug_par_each_eval = TRUE, debug_par_full = TRUE,
    debug_par_max_eval = Inf, tol_loglik = 1e-07, tol_par = 0.001, Q = NULL,
    simple = NULL, nquadgh_u = nquadgh, nquadgh_zeta = nquadgh, nquadgh_re = nquadgh,
    rho_pattern = NULL, two_tier = NULL)
}
\arguments{
\item{data}{A data frame with one row per subject and occasion. Use integer categories starting at 0 or 1, or ordered factors. Use NA for missing item responses.}
\item{namvar}{Names of the item columns.}
\item{occ, ID}{Names of the occasion and subject columns. The occasion labels determine time order.}
\item{identification, marker_item}{Identification scheme and, for marker or mixed_marker, one active marker item per factor.}
\item{intercept_time_free}{Items whose intercepts can vary after baseline, or an item by occasion indicator matrix.}
\item{intercept_covariates, slope_covariates}{Names of subject-level covariates that predict the growth intercept and slope.}
\item{rho_fixed}{Optional fixed score in the original three-occasion rho model. A vector can specify one score per factor. Use \code{rho_pattern} to set scores at arbitrary occasions.}
\item{rho_pattern}{A numeric vector in sorted occasion order or a matrix with one row per occasion and one column per factor. The vector is repeated across factors. Each factor must have a fixed score of 0 and a fixed score of 1. Each \code{NA} defines a separate estimated score. Requires the simple growth design and selects \code{growth_shape = "rho"} when \code{growth_shape} is omitted.}
\item{methInteg}{Integration rule. The general fitter supports GH, QMC, AGH, AGH_FULL_EQUIV, AGH_RE_NESTED and AGH_RE_FULL_EQUIV.}
\item{two_tier}{Use separate rules for the latent trajectory and item intercepts. Set to \code{FALSE} for standard GH or AGH. QMC uses joint integration and requires \code{FALSE}. The default \code{NULL} selects joint QMC and two-tier GH or AGH. \code{AGH_FULL_EQUIV} with \code{FALSE} selects standard AGH. The growth coefficient AGH rules use the two-tier layout.}
\item{re_item, theta_time}{Include Gaussian item random intercepts and allow item residual variances to vary by occasion.}
\item{monotone_tau}{Parameterize threshold increments to be positive.}
\item{loadmask, factor_names, Q, simple}{A 0/1 matrix of active loadings and optional factor names. Q and simple can specify the number of factors and one loading per item when loadmask is omitted.}
\item{threshold_identification}{Fix the first two or central two thresholds to zero and one. Binary items use one fixed threshold and unit residual variance.}
\item{growth_shape}{Use an additive change to the selected time score, or specify fixed and estimated scores with \code{rho_pattern}. The original three-occasion rho specification is also available.}
\item{nquadgh, nquadgh_u, nquadgh_zeta, nquadgh_re}{Gauss-Hermite nodes per dimension from \pkg{fastGHQuad}. Separate orders are available for growth effects, occasion disturbances and item intercepts in the two-tier rules.}
\item{nquadqmc}{Number of Sobol points for QMC.}
\item{update_gh}{Use zero for a refreshed adaptive grid at each evaluation. A positive value gives the initial number of optimizer iterations per fixed grid.}
\item{analytic_score}{Use the analytic score in supported original three-occasion rho models with one factor, two random coefficients and a fixed adaptive grid. Models with \code{rho_pattern} use numerical derivatives.}
\item{delta_free}{Estimate the adjustment to the masked time score in the simple offset model.}
\item{simple_growth, same_design, subject_specific_design}{Use the simple intercept and slope design, share designs across factors, and allow designs to vary by subject, respectively.}
\item{time_scores}{Scores for the simple design, possibly varying by subject or factor. An occasion vector or occasion by factor matrix with \code{NA} entries can be used in place of \code{rho_pattern}.}
\item{fixed_covariates, time_interactions}{Covariates and their optional interactions with time in the simple growth design.}
\item{X, Z}{Fixed and random designs for simple_growth = FALSE. A common design has one row per occasion. Subject-specific designs may use arrays indexed by subject, occasion, factor and column. Subject and occasion dimensions use sorted labels.}
\item{delta_mask}{An indicator selecting the scores changed by delta or rho in the simple design.}
\item{fixed_names, random_names}{Names of the fixed and random design columns.}
\item{maxiter, start, parscale}{Iteration budget, optional starting vector and parameter scaling. Use zero iterations to evaluate the likelihood at the starting values.}
\item{compute_hessian, hessian_method, hess_eps}{Whether to compute the observed information, derivative method, and step for the analytic-score Hessian. Joint integration uses numerical derivatives.}
\item{profile, profile_print, verbose}{Collect timings, print timings, and display optimization progress.}
\item{debug_par, debug_par_each_eval, debug_par_full, debug_par_max_eval}{Set \code{debug_par = TRUE} to print parameter diagnostics. The other arguments control output at each evaluation, full parameter vectors, and the evaluation limit.}
\item{tol_loglik, tol_par}{Relative convergence tolerances for the likelihood and coefficients.}
}
\details{
The fixed and random design matrices may have different numbers of columns. Each column of Z gives one Gaussian random growth coefficient per factor. Occasion disturbances are independent over occasions and factors.
The simple design adds an intercept and a slope for each factor. The original rho specification uses occasions coded 0, 1, and 2. For five occasions, \code{rho_pattern = c(0, 1, NA, NA, NA)} estimates the scores at occasions three to five. The pattern \code{c(0, NA, NA, NA, 1)} estimates the scores at occasions two to four. The fixed zero anchors the intercept and the fixed one sets the slope scale. The time basis can be specified separately by factor. The custom design uses the supplied X and Z matrices.

Effects identification sets the mean of active loadings to one and the sum of their item intercepts to zero within each factor. Baseline standardization fixes the growth intercept mean to zero and the variance at its zero score to one at reference covariate values. This requires the simple design. Marker identification fixes the marker loading to one and its baseline item intercept to zero. Mixed marker identification fixes the marker loading to one and the growth intercept mean to zero. The time anchors identify the time basis. The full model requires sufficient observations and independent design columns.

FIML uses the observed item categories for each subject. Let \eqn{T} be the number of occasions, \eqn{Q} the number of factors and \eqn{K} the number of included item intercepts. Standard GH and AGH use \eqn{m^{TQ + K}} joint nodes at order \eqn{m}. QMC draws the specified number of Sobol points in \eqn{TQ + K} dimensions. Two-tier GH and trajectory AGH integrate the trajectory in \eqn{TQ} dimensions. Each included item intercept contributes a separate one-dimensional rule. The full-equivalent adaptive rule uses the joint mode and Hessian. The growth coefficient rules integrate \eqn{Qp} growth coordinates and \eqn{TQ} occasion disturbances, with \eqn{p} random coefficients per factor. See \code{\link{integration_info}} for dimensions and node counts.
}
\value{A fitted object with the log likelihood, convergence code, parameter estimates and decoded measurement and growth parameters. The \code{BIC} field uses the total number of subjects. The \code{decode} field contains \code{beta}, \code{Sigma_u} and \code{zeta_var}. With a time-score pattern, \code{rho} contains the occasion, factor, score, estimated indicator and standard error. The standard error of an estimated score is available when the information matrix can be inverted.}
\examples{
dat <- expand.grid(ID = 1:12, occ = 0:4)
dat$a <- rep(1:4, length.out = nrow(dat))
dat$b <- rep(c(2L, 4L, 1L, 3L, NA_integer_, 2L), length.out = nrow(dat))
B <- cbind(1, 0:4, (0:4)^2)
fit <- growth_tt_fit_flexible(dat, c("a", "b"), simple_growth = FALSE,
                              X = B, Z = B, nquadgh = 2L, maxiter = 0L, debug_par = FALSE)
dim(fit$decode$Sigma_u)

rho_fit <- growth_tt_fit(dat, c("a", "b"), rho_pattern = c(0, NA, NA, NA, 1),
                         nquadgh = 2L, maxiter = 0L, debug_par = FALSE)
rho_fit$rho
integration_info(rho_fit)
}
