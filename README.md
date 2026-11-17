# LongitudinalLVM

The R package `longProbitLVM` fits longitudinal ordinal probit models by full information maximum likelihood (FIML). Supply one row per subject and occasion, with `NA` for missing responses. Items can be integer coded from 0 or 1, or stored as ordered factors.

## Installation

Install from GitHub with `remotes::install_github("niccolocao/LongitudinalLVM")`. From a local copy, install `fastGHQuad` and run `R CMD INSTALL .` in the source directory. The package compiles C and Fortran code and links to R's BLAS and LAPACK libraries. Install `spacefillr` to use QMC.

## General growth model

Use `growth_tt_fit_flexible()` with `simple_growth = FALSE` and design matrices `X` and `Z`. Their columns define the fixed and random coefficients. Common designs have one row per occasion, in sorted occasion order. The function also accepts subject-specific designs and separate designs by factor. The following example fits a quadratic trajectory.

```r
library(longProbitLVM)

s <- sort(unique(dat$measurecount))
s <- s - s[1L]
B <- cbind(intercept = 1, linear = s, quadratic = s^2)

fit <- growth_tt_fit_flexible(
  dat, item_names, occ = "measurecount",
  simple_growth = FALSE, X = B, Z = B, delta_free = FALSE,
  methInteg = "AGH_FULL_EQUIV", nquadgh = 3L
)
fit$decode$beta
fit$decode$Sigma_u
```

In this example, the population means and the individual trajectories include a quadratic term. Use a different `Z` to choose which terms vary between subjects. For example, use `Z = B[, 1:2]` for a fixed quadratic and random intercept and linear slope.

For factor $q$ and subject $i$, the growth equation is

$$
\eta_{itq}=\mathbf x_{itq}^{\mathsf T}\boldsymbol\beta_q+
\mathbf z_{itq}^{\mathsf T}\mathbf u_{iq}+\zeta_{itq}.
$$

Stacking all factors and occasions gives
$\boldsymbol\eta_i\sim N(X_i\boldsymbol\beta,\,Z_i\Sigma_u Z_i^{\mathsf T}+D_\zeta)$.
The covariance $\Sigma_u$ is free across growth coefficients and factors. The occasion disturbances have a diagonal covariance $D_\zeta$. With a common design, three occasions, and three random coefficients per factor, a full random coefficient covariance and separate occasion variances are confounded.

| Function or option | Design |
| --- | --- |
| `simple_growth = TRUE` | Random intercept and slope, with default scores $0,1,\ldots,T-1$. Named covariates can predict the intercept and slope. |
| `growth_shape = "offset"` | The default mask adds $\delta_q$ to the final score for factor $q$. Use `delta_free = FALSE` to keep the specified scores fixed. |
| `growth_shape = "rho"` | Use `rho_pattern` for any placement of fixed and free scores. In the original three-occasion model, set the third score with `rho_fixed`. |

### Fixed and estimated time scores

Use one entry per sorted occasion. Supply a vector shared by the factors or a matrix with one column per factor. Each column must contain a fixed $0$ and a fixed $1$. All other numeric entries are fixed. Each NA is a separate parameter.

$$
\eta_{itq}=(\beta_{0q}+u_{0iq})+(\beta_{1q}+u_{1iq})s_{tq}+\zeta_{itq}.
$$

The zero score defines the intercept. The difference between the scores zero and one defines the unit of the slope. Values at free positions are unconstrained.

```r
fit <- growth_tt_fit(
  dat, item_names, occ = "measurecount",
  rho_pattern = c(0, NA, NA, NA, 1),
  methInteg = "AGH_FULL_EQUIV", nquadgh = 3L
)
fit$rho
integration_info(fit)
```

With `c(0, 1, NA, NA, NA)`, occasions 1 and 2 fix the time scale. With `c(0, NA, NA, NA, 1)`, occasions 1 and 5 fix it. A free row in `fit$rho` includes a standard error when `compute_hessian = TRUE` yields an invertible positive definite information matrix. Missing item responses are entered as NA in `data`.

### Growth identification

| `identification` | Factor scale and location |
| --- | --- |
| `"effects"` | Active loadings average one and the corresponding item intercepts sum to zero for each factor. |
| `"baseline_std"` | The growth intercept mean is zero. The latent variance at score zero is one at reference covariate values. This requires the simple design. |
| `"marker"` | The marker loading is one and its baseline item intercept is zero. The growth intercept mean is estimated. |
| `"mixed_marker"` | The marker loading is one and the growth intercept mean is zero. |

For growth models, `threshold_identification` fixes the first two (`"first_two"`) or central two (`"central_01"`) thresholds at zero and one. Binary items use a zero threshold and unit residual variance. Growth loadings and thresholds are common across occasions. Use `intercept_time_free` for selected time-specific item intercepts and `theta_time` for time-specific residual variances. The two time anchors fix the location and scale of the time basis. Estimability of the full model also depends on the measurement constraints, observed occasion overlap and the rank of the fixed and random designs.

## Measurement and CFA

For item $k$ at time $t$, let
$m_{ikt}=\nu_{kt}+\lambda_{kt}^{\mathsf T}\eta_{it}+a_{ik}$, where $a_{ik}\sim N(0,\psi_k)$ is an optional subject by item random intercept. The category probability is

$$
p_{ikt}(c)=\Phi\!\left(\frac{\tau_{kt,c}-m_{ikt}}{\sqrt{\theta_{kt}}}\right)
-\Phi\!\left(\frac{\tau_{kt,c-1}-m_{ikt}}{\sqrt{\theta_{kt}}}\right).
$$

The residuals are conditionally independent. Their variances $\theta_{kt}$ can be common across occasions or time specific. The shared $a_{ik}$ induces dependence between repeated responses to the same item.

- `cfa_long_fit_md()` estimates free means and covariance for the stacked latent factors subject to identification constraints. It supports covariates in the latent mean, cross loadings and selected measurement invariance constraints. Each item requires at least four categories. Fix one loading per factor at one and use `threshold_id` to anchor thresholds 1 and 2 or 2 and 3. With `mean_ident = "latent_baseline"`, anchor the first occasion's factor means at zero. With `mean_ident = "marker_intercept"` in full-equivalent two-tier AGH, anchor the marker intercepts at zero. The order in `meas_invar` is thresholds, intercepts, loadings and residual variances.
- `cfa_long_fit_wu2016()` fits a selected invariance stage. The stages successively constrain thresholds, loadings, intercepts and residual variances. Each item requires one active loading and at least three categories. `measurement_invariance_wu2016()` fits and compares the stages. [Wu and Estabrook (2016)](https://doi.org/10.1007/s11336-016-9506-0)

## Observed data likelihood and numerical integration

For a subject $i$, let $O_{ik}$ be the observed occasions for item $k$. The contribution is

$$
L_i=E_{\eta_i}\!\left[\prod_k E_{a_{ik}}\!\left[
\prod_{t\in O_{ik}}p_{ikt}(y_{ikt})
\right]\right],\qquad
\ell=\sum_i\log L_i.
$$

An empty product is one. The missingness mechanism is ignorable for likelihood inference under MAR and distinct parameters. [Rubin (1976)](https://doi.org/10.1093/biomet/63.3.581)

For joint integration on standard normal coordinates $z$, write
$g_i(z)=\prod_{(k,t)\in O_i} p_{ikt}(y_{ikt}\mid z)$.

Standard GH uses $\widehat L_i=\sum_j w_jg_i(z_j)$. The package generates the nodes and weights with `fastGHQuad::gaussHermiteData()` and scales the rule for a standard normal integral. QMC uses $\widehat L_i=M^{-1}\sum_jg_i(\Phi^{-1}(v_j))$ with Owen-scrambled Sobol points $v_j$ and the inverse CDF applied coordinate-wise.

For joint AGH, let $z_{ij}^*=\hat z_i+C_i z_j$, with $C_iC_i^{\mathsf T}=H_i^{-1}$, where $H_i$ is the Hessian of the negative log integrand including the Gaussian density. Then

$$
\widehat L_i=|\det C_i|\sum_j w_jg_i(z_{ij}^*)
\exp\!\left(\frac{\|z_j\|^2-\|z_{ij}^*\|^2}{2}\right).
$$

In two-tier integration, each outer node fixes the trajectory. Apply each item rule to its own product of observed response probabilities, then multiply the summed item contributions. For GH this gives

$$
\widehat L_i=\sum_j w_j^\eta \prod_k\left[
  \sum_h w_h^a \prod_{t\in O_{ik}}p_{ikt}(y_{ikt}\mid\eta_j,a_{kh})
\right].
$$

In AGH, the node-specific rule adapts each item integral conditional on the outer node. The full-equivalent rule uses a block factor of the joint inverse Hessian with the same item separation in the sum. [Pinheiro and Bates (1995)](https://doi.org/10.1080/10618600.1995.10474663)

### Integration dimensions

Let $T$ count occasions, $Q$ factors, $K$ included item random intercepts and $p$ random growth coefficients per factor. Write $d=TQ$ and $b=Qp$. When `re_item = FALSE`, $K=0$ and item rules are omitted. An $r$-dimensional GH or AGH rule of order $m$ uses $m^r$ nodes. Write $m_u$, $m_\zeta$ and $m_a$ for separate orders in the growth coefficient representation.

| Algorithm | Integration blocks | Largest mode and Hessian | Configured nodes |
| --- | --- | --- | --- |
| Standard GH | Joint $d+K$ | None | $m^{d+K}$ |
| Standard AGH | Joint $d+K$ | $d+K$ | $m^{d+K}$ |
| Two-tier GH | Outer $d$, then $K$ rules of dimension 1 | None | $m^d$ outer, $m$ per item at each outer node |
| Two-tier AGH, node-specific | Outer $d$, then $K$ rules of dimension 1 | $d$ | $m^d$ outer, $m$ per item at each outer node |
| Two-tier AGH, full-equivalent | Outer $d$, then $K$ rules of dimension 1 | $d+K$ | $m^d$ outer, $m$ per item at each outer node |
| `AGH_RE_NESTED` | Outer $b$, then $d$ occasion disturbances and $K$ rules of dimension 1 | $b$ | $m_u^b$ adaptive nodes, $m_\zeta^d$ fixed nodes per adaptive node, $m_a$ per item |
| `AGH_RE_FULL_EQUIV` | Outer $b+d$, then $K$ rules of dimension 1 | $b+d+K$ | $m_u^b m_\zeta^d$ outer, $m_a$ per item at each outer node |
| QMC | Joint $d+K$ | None | $M$ Sobol points in $d+K$ dimensions |

In `growth_tt_fit()` and `cfa_long_fit_md()`, select `two_tier = FALSE` with `methInteg = "GH"` or `"AGH"` for standard integration. QMC uses one joint rule with `two_tier = FALSE`. Use `two_tier = TRUE` for two-tier GH or AGH. In CFA, `two_tier_type` selects `"node_specific"` or `"full_equiv"` for two-tier AGH. In growth, `"AGH"` and `"AGH_FULL_EQUIV"` select those two rules when `two_tier = TRUE`. The growth coefficient methods let you set `nquadgh_u`, `nquadgh_zeta` and `nquadgh_re` separately. Wu CFA uses a trajectory rule with $K=0$. `integration_info(fit)` reports the dimensions and counts, including the dimension of the adaptive mode.

The trajectory and growth coefficient representations use the same model likelihood. Their finite quadrature rules can give different numerical values. `growth_tt_individual_derivatives()` evaluates subject likelihood contributions and finite difference derivatives for GH, QMC and trajectory-based AGH. The analytic score is available for supported single factor, two random coefficient models with the original three-occasion `growth_shape = "rho"` specification. The fitter uses numerical derivatives for `rho_pattern`.

## Reproducing the original growth fit

For the original three-occasion fit with estimated scores $(0,1,\rho)$, use `growth_shape = "rho"`, `threshold_identification = "central_01"`, `methInteg = "AGH_FULL_EQUIV"` and `two_tier = TRUE`. Match the original quadrature order and the identification and covariate options. The quadrature uses the original node scaling, tensor order and weight products. Set `debug_par = FALSE` to suppress parameter diagnostics. The fitted object includes the data report, profiling summary and the original BIC based on the total number of subjects.
