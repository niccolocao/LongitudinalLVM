# Cross-sectional ordinal probit CFA (experimental standalone module)

This module adds cross-sectional CFA without modifying any longitudinal code,
package metadata or namespace. It is installed as a supplementary R script.
It can also be sourced directly from this directory.

## Installation and loading

    remotes::install_github("niccolocao/LongitudinalLVM",
                            ref = "cross-sectional-cfa")
    source(system.file("cross-sectional", "cfa_fit.R",
                       package = "LongitudinalLVM"))

No change to NAMESPACE or DESCRIPTION is required. The function is not
an exported package function; call cfa_fit() after sourcing the file.

## Verbose progress

Set verbose = TRUE to display the method, sample dimensions, number of
quadrature nodes, starting likelihood, BFGS progress every ten optimizer
iterations (REPORT = 10), convergence code and post-estimation
calculations. Initial, final and AGH outer-cycle messages are also shown. This is available for GH, AGH and QMC.
The default verbose = FALSE preserves silent fitting.

    fit <- cfa_fit(dat, Q = 1, methInteg = "GH",
                   nquadgh = 50, verbose = TRUE)

    fit_agh <- cfa_fit(dat, Q = 1, methInteg = "AGH",
                       nquadgh = 20, verbose = TRUE)

Objective evaluations include numerical-gradient evaluations and are not
the same as optimizer iterations. AGH progress is reported for the outer
optimization, not every individual posterior mode search. Verbose output
does not alter the likelihood or optimization settings.

## Model

For subject i and item k,

    Y*_{ik} = nu_k + lambda_k' eta_i + epsilon_{ik}
    eta_i ~ N_Q(mu, Phi)
    epsilon_{ik} ~ N(0, 1)

Residuals are mutually independent and independent of the factors.
The observed category is determined by the ordered thresholds.
The first finite threshold is zero for every item. All residual
variances are one. Additional finite thresholds are freely estimated
in increasing order. For each factor, one pure marker loading is fixed
to one and the corresponding item intercept is fixed to zero.
All other active loadings, nonmarker intercepts, factor means and the
positive definite factor covariance matrix are free.

This is a marker-identified model, not an unrestricted exploratory factor
model. Each factor requires a pure marker item. The loading pattern must
also identify the remaining parameters; a pure marker alone does not ensure
identification. An information-rank warning is reported when detected.

## Single factor, five binary items

    set.seed(123)
    N <- 2000
    eta <- rnorm(N, 0.3, sqrt(1.2))
    lambda <- c(1, .8, 1.2, .9, 1.1)
    nu <- c(0, .3, -.3, .2, -.2)
    ystar <- sweep(outer(eta, lambda), 2, nu, "+") +
      matrix(rnorm(N * 5), N, 5)
    dat <- as.data.frame(1L * (ystar > 0))
    names(dat) <- paste0("y", 1:5)

    fit <- cfa_fit(dat, Q = 1, marker_item = "y1",
                   methInteg = "GH", nquadgh = 50,
                   compute_hessian = TRUE)

    fit$coefficients
    fit$decode
    fit$convergence
    fit$loglik
    min(N * fit$probabilities)
    colSums(fit$scores)

The ten coefficients are ordered as the four free loadings, four free
intercepts, factor mean and factor variance. Natural-scale names are
lambda[y2,F1], ..., nu[y2], ..., mu[F1], phi[F1,F1].
For ordinal items, additional threshold parameters follow the covariance
parameters.

## Adaptive Gauss-Hermite quadrature (AGH)

    fit_agh <- cfa_fit(dat, Q = 1, marker_item = "y1",
                       methInteg = "AGH", nquadgh = 12,
                       compute_hessian = TRUE)

    fit_agh$integration
    fit_agh$convergence
    fit_agh$loglik

AGH locates the posterior mode of the standardized latent factor vector
separately for every individual, including subjects with identical response
patterns. No pattern-level caching is used during AGH estimation or observed
score and Hessian evaluation. The negative posterior
log-density Hessian determines the covariance of the adapted normal
quadrature rule. The change of variables includes its Jacobian and the
normal density ratio. Missing items are omitted from the conditional
likelihood; entirely missing subjects have probability one and require no mode search.

By default, update_gh = 10 selects the fast outer-inner AGH algorithm
based on the longitudinal functions. At the beginning, an adaptive grid
is constructed for each individual in standardized latent coordinates.
The inner BFGS optimization keeps these nodes and importance corrections
fixed for a block of at most ten iterations, including numerical gradient
evaluations. The outer step refreshes all individual grids, warm-starting
the mode search at the previous individual mode. Convergence requires
relative changes of both the refreshed likelihood and parameters below
tol_loglik (1e-7) and tol_par (1e-3). Blocks have constant length, unlike
the increasing inner budgets in the original longitudinal functions.

Set update_gh = 0 for full AGH, which refreshes the grid at every objective
evaluation. Final observed-data likelihood, numerical scores, Jacobian
and optional Hessian use refreshed adaptive quadrature. The numerical
derivatives can still be expensive. The implementation is in R with
vectorized fixed-grid likelihood evaluation, not the longitudinal
compiled Fortran routines.

The mode algorithm uses BFGS with an analytic posterior gradient and
Newton refinement. The controls agh_maxiter and agh_tol set the
mode-search iteration budget and gradient tolerance. An unsuccessful
adaptation is not silently replaced with standard GH.

    fit_agh <- cfa_fit(dat, Q = 1, marker_item = "y1",
                       methInteg = "AGH", nquadgh = 20,
                       update_gh = 10, verbose = TRUE)

    fit_agh$agh_trace
    fit_agh$integration$update_gh

The trace records inner iteration budgets, fixed and refreshed
likelihoods, convergence measures and actual individual mode searches.
The inner mode-search count must be zero; the refresh count equals the
number of individuals with observed responses. Early BFGS convergence
can end an inner block before its maximum iteration budget.

For a Q-dimensional latent vector, AGH uses nquadgh^Q nodes per
individual with observed responses. There is no fixed factor limit, but max_nodes applies
to both GH and AGH. Adaptation can reduce the required quadrature order,
although the mode searches introduce overhead. In particular, repeated
response patterns are not aggregated for AGH, so fitting and numerical
derivatives can be much slower than GH for large samples.

    fit_agh$integration$adaptation
    fit_agh$integration$final_adaptations
    fit_agh$individual_loglik

The final_adaptations counter records the number of mode searches in a
single final observed-data likelihood evaluation. It equals the number
of subjects with at least one observed item. Complete-pattern probability
and Jacobian outputs are still evaluated separately on their requested
pattern grid, not on the observed individuals.

For validation, compare with a high-order GH reference at identical
parameters before comparing optimized fits:

    ref <- cfa_fit(dat, start = fit_agh$coefficients, maxiter = 0,
                   methInteg = "GH", nquadgh = 50)
    max(abs(ref$probabilities - fit_agh$probabilities))
    ref$loglik - fit_agh$loglik
    max(abs(ref$scores - fit_agh$scores))

Unlike a common fixed quadrature grid, pattern-specific AGH rules
need not produce complete-pattern probabilities summing exactly to one
at finite quadrature order. The probabilities are not renormalized;
increase nquadgh and verify normalization, derivatives and likelihood
stability before using them for goodness-of-fit inference.

## Multiple factors

    fit2 <- cfa_fit(dat2, Q = 2,
                    simple = c(1, 1, 1, 2, 2, 2),
                    marker_item = c("y1", "y4"),
                    nquadgh = 25)

Use simple for simple structure, or specify a K by Q binary loadmask
to allow cross-loadings. For example:

    mask <- rbind(c(1, 0), c(1, 0), c(1, 1),
                  c(0, 1), c(0, 1), c(1, 1))
    fit2 <- cfa_fit(dat2, loadmask = mask,
                    marker_item = c("y1", "y4"),
                    nquadgh = 25)

The covariance matrix is unrestricted positive definite, using a
log-Cholesky parameterization internally. There is no fixed maximum
number of factors, but tensor Gauss-Hermite quadrature grows as
nquadgh^Q. The max_nodes control prevents an excessive node count.
For higher dimensions, fixed quasi-Monte Carlo integration is available:

    set.seed(123)
    fit_qmc <- cfa_fit(dat2, simple = c(1, 1, 1, 2, 2, 2),
                       methInteg = "QMC", nquadqmc = 32768)

QMC requires the optional spacefillr package. Set a seed for reproducibility
and check stability by increasing the number of nodes. Numerical integration
error must be assessed before using either integration rule for inference.

## Outputs

| Field | Definition |
| --- | --- |
| coefficients | Free parameters on the natural scale |
| decode | Loading matrix, intercepts, latent means/covariance, thresholds and residual variances |
| loglik | Maximized observed-data marginal log likelihood (if converged) |
| individual_loglik | Observed-data log-likelihood contributions in input row order |
| convergence | optim convergence code; zero indicates optimizer convergence |
| patterns | Complete response patterns, zero-based categories |
| probabilities | Model-implied probabilities for complete patterns |
| jacobian | Derivatives of complete pattern probabilities with respect to natural parameters |
| pattern_scores | Derivatives of log probabilities for complete patterns |
| scores | Individual observed-data likelihood scores, in input row order |
| information | Expected full-data Fisher information per subject |
| empirical_information | Outer product of observed-data scores divided by N |
| hessian | Optional total observed negative log-likelihood Hessian, natural scale |
| observed_patterns | Distinct observed patterns, including missingness |
| observed_counts | Frequencies of the observed patterns |
| probability_function | Evaluate probabilities at supplied natural parameters and patterns |
| jacobian_function | Evaluate probability derivatives at supplied natural parameters and patterns |
| integration | Integration method, dimension, base node count, adaptation indicator and final adaptation count for AGH |

The probability and derivative functions use the same integration method
as estimation. GH and QMC reuse a fixed numerical rule; AGH recomputes
the individual-specific adaptation at each observed-data evaluation. Missing responses are marginalized under
the observed-data likelihood and usual ignorability assumptions. Numeric
items must have consecutive categories starting at 0 or 1, or be ordered
factors. All declared categories must be observed for estimation. Patterns
in the outputs and evaluation functions always use zero-based codes; the
levels field provides the mapping to original category labels.

Full pattern enumeration is optional and limited by max_patterns.
Set full_patterns = FALSE for large item sets. Individual scores are
still available. Expected complete-data information is not the
observed-data information when responses are missing.

For simulations, start accepts a natural-scale numeric vector in the
coefficient order, or a vector named exactly as coefficients. Set maxiter
to zero to evaluate at start without optimizing. The resulting convergence
code is NA. Derivatives use adaptive-step central finite differences.
Hessian computation is optional and may fail at parameter boundaries.

## Goodness-of-fit example (complete binary data)

    id <- match(apply(as.matrix(dat), 1, paste, collapse = ":"),
                apply(fit$patterns, 1, paste, collapse = ":"))
    n <- tabulate(id, nrow(fit$patterns))
    expected <- sum(n) * fit$probabilities
    X2 <- sum((n - expected)^2 / expected)
    pos <- n > 0
    G2 <- 2 * sum(n[pos] * log(n[pos] / expected[pos]))
    df <- length(n) - 1 - fit$npar
    pchisq(c(X2, G2), df, lower.tail = FALSE)

The usual chi-square calibration requires an identified, interior,
correctly specified model and adequate expected frequencies. This
module supplies estimation and derivatives, not the sign-flipping tests.

## Validation

See the new tests/cfa-cross-sectional.R script in the repository.
The module is experimental until the R tests and numerical comparisons
have been run. No longitudinal functions have been changed.
