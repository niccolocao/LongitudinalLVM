# LongitudinalLVM

`LongitudinalLVM` fits longitudinal latent variable models for ordinal responses using a probit link. It supports longitudinal confirmatory factor analysis (CFA) and second-order latent growth models. Estimation uses full-information maximum likelihood (FIML), with several numerical integration methods.

## Installation

```r
remotes::install_github("niccolocao/LongitudinalLVM")
```

The package compiles C and Fortran code and links to BLAS and LAPACK. Install `fastGHQuad` for Gauss-Hermite quadrature and `spacefillr` for quasi-Monte Carlo integration.

## Model specification

### Measurement model

Let $Y_{ijk}$ be the ordinal response of subject $i$ to item $k$ at occasion $j$, with categories $c=0,\ldots,C_k-1$. Following the underlying variable approach,

$$
Y_{ijk}=c \iff \tau_{jk,c}\leq Y^*_{ijk}<\tau_{jk,c+1},
$$

where $\tau_{jk,0}=-\infty$ and $\tau_{jk,C_k}=+\infty$. The measurement model is

$$
Y^*_{ijk}=\nu_{jk}+\lambda_{jk}\eta_{ij}+b_{ik}+\epsilon_{ijk},
\qquad \epsilon_{ijk}\sim N(0,\sigma^2_{\epsilon,jk}),
$$

where $\eta_{ij}$ is the occasion-specific latent factor and $b_{ik}\sim N(0,d_k^2)$ is an optional item-specific random intercept shared across occasions. The residuals are conditionally independent. For multiple factors, $\lambda_{jk}\eta_{ij}$ is replaced by $\boldsymbol\lambda_{jk}^{\top}\boldsymbol\eta_{ij}$.

The response probabilities are obtained from the normal CDF $\Phi$:

$$
\Pr(Y_{ijk}=c\mid\eta_{ij},b_{ik})=
\Phi\!\left(\frac{\tau_{jk,c+1}-m_{ijk}}{\sigma_{\epsilon,jk}}\right)
-\Phi\!\left(\frac{\tau_{jk,c}-m_{ijk}}{\sigma_{\epsilon,jk}}\right),
$$

with $m_{ijk}=\nu_{jk}+\lambda_{jk}\eta_{ij}+b_{ik}$ for a single factor.

Longitudinal measurement invariance can be imposed on thresholds, loadings, and item intercepts. Residual variances may be occasion-specific. The options `intercept_time_free` and `theta_time` control selected time-specific intercepts and residual variances.

### Second-order latent growth

For a single latent factor, the structural model is

$$
\eta_{ij}=\upsilon_{i0}+z_j\upsilon_{i1}+\zeta_{ij},
$$

where $\upsilon_{i0}$ and $\upsilon_{i1}$ are the latent intercept and slope, and $\zeta_{ij}\sim N(0,\omega_j)$ is an occasion-specific disturbance. Covariates can predict the growth factors, whose covariance is freely estimated.

With three occasions, the linear model uses $\mathbf z=(0,1,2)^\top$; the latent-basis model uses $\mathbf z=(0,1,\rho)^\top$, with $\rho$ estimated. The package also supports multiple factors and more general fixed and random growth designs.

| Option | Growth specification |
| --- | --- |
| `simple_growth = TRUE` | Random intercept and slope; default scores $0,1,\ldots,n-1$ |
| `growth_shape = "offset"` | Estimated offset to the final score; `delta_free = FALSE` fixes it |
| `growth_shape = "rho"` | Fixed and estimated scores specified by `rho_pattern`; `rho_fixed` fixes the third score in the original three-occasion model |

In `rho_pattern`, each `NA` is an estimated score. Each factor must have fixed scores 0 and 1 to identify the time scale.

### Identification

For growth models, `identification` accepts `"effects"`, `"baseline_std"`, `"marker"`, and `"mixed_marker"`. The `"marker"` specification fixes a marker loading to one and its baseline item intercept to zero. `"baseline_std"` is available for the simple growth design.

`threshold_identification` fixes either the first two (`"first_two"`) or central two (`"central_01"`) finite thresholds at 0 and 1. Binary items use a zero threshold and unit residual variance.

## Model fitting

- `growth_tt_fit()` fits second-order latent growth models with covariates and flexible time scores.
- `cfa_long_fit_md()` fits longitudinal CFA models with latent means, covariance parameters, covariates, cross-loadings, and selected invariance constraints. It requires at least four categories per item.
- `cfa_long_fit_wu2016()` fits a longitudinal CFA invariance stage; `measurement_invariance_wu2016()` fits and compares stages following [Wu and Estabrook (2016)](https://doi.org/10.1007/s11336-016-9506-0). These methods require at least three categories per item.

### FIML and numerical integration

Under missing at random (MAR) and the usual ignorability conditions, FIML uses all observed responses. Let $O_{ik}$ be the observed occasions for subject $i$ and item $k$. The observed-data likelihood separates item-specific random effects into one-dimensional integrals conditional on the latent factors:

$$
L_i=\int\prod_{k=1}^K\left[
\int\prod_{j\in O_{ik}}
\Pr(Y_{ijk}=y_{ijk}\mid\eta_{ij},b_{ik})
\,h_b(b_{ik})\,db_{ik}\right]
h_\eta(\boldsymbol\eta_i\mid\mathbf w_i)\,d\boldsymbol\eta_i.
$$

Available integration methods include standard Gauss-Hermite (GH), adaptive Gauss-Hermite (AGH), quasi-Monte Carlo (QMC), and two-tier GH/AGH. Two-tier methods integrate over latent factors in an outer rule and over each item-specific random effect in a separate inner rule. Two-tier AGH supports node-specific adaptation and full-equivalent adaptation based on the joint posterior approximation.

Use `methInteg` to select the method and `two_tier` to switch between joint and two-tier integration. Growth coefficient methods also support `AGH_RE_NESTED` and `AGH_RE_FULL_EQUIV`, with separate quadrature orders set through `nquadgh_u`, `nquadgh_zeta`, and `nquadgh_re`. `integration_info(fit)` reports the integration dimensions and node counts.

## Example: estimated time scores

```r
fit <- growth_tt_fit(
  dat, item_names, occ = "measurecount",
  rho_pattern = c(0, NA, NA, NA, 1),
  methInteg = "AGH_FULL_EQUIV", nquadgh = 3L
)

fit$rho
integration_info(fit)
```

Missing responses are supplied as `NA`. For the original three-occasion latent-basis model with scores $(0,1,\rho)$, use `growth_shape = "rho"`, `threshold_identification = "central_01"`, `methInteg = "AGH_FULL_EQUIV"`, and `two_tier = TRUE`, together with the matching identification, covariates, and quadrature order.
