\name{growth_tt_individual_derivatives}
\alias{growth_tt_individual_derivatives}
\title{Calculate subject likelihood contributions and derivatives}
\description{Calculate log likelihood contributions and finite difference derivatives for GH, QMC and adaptive quadrature over the latent trajectory.}
\usage{
growth_tt_individual_derivatives(fit, data = NULL, namvar = fit$c_inputs$namvar, occ = "occ",
    ID = "ID", nquadgh = NULL, nquadqmc = 2000, eps_score = 1e-05, eps_hess = 1e-04,
    par_index = NULL, compute_score = TRUE, compute_hessian = TRUE, agh_fixed_grid = TRUE,
    progress = FALSE, attach_to_fit = FALSE)
}
\arguments{
\item{fit}{A fitted ordinal growth model.}
\item{data, namvar, occ, ID}{Optional data and column names. Fits with stored inputs use those inputs to evaluate the likelihood.}
\item{nquadgh, nquadqmc}{Quadrature order or number of QMC points when reconstructing the inputs of an earlier fit.}
\item{eps_score, eps_hess}{Relative finite difference step sizes.}
\item{par_index}{Optional indices of the coefficients to differentiate with respect to.}
\item{compute_score, compute_hessian}{Compute subject scores and Hessians.}
\item{agh_fixed_grid}{Keep the fitted adaptive grid fixed when taking derivatives. Required for the adaptive rules.}
\item{progress}{Print derivative evaluation progress.}
\item{attach_to_fit}{Attach results to the fit and return it.}
}
\details{
Subjects with no observed items contribute zero to the log likelihood.
}
\value{A list containing subject log likelihoods, optional derivatives, a sandwich estimate when available and numerical checks, or the fit with these results attached.}
