\name{integration_info}
\alias{integration_info}
\title{Report the dimensions of numerical integration}
\description{Report the configured integration blocks, mode dimensions and node counts for a fitted model.}
\usage{
integration_info(fit)
}
\arguments{
\item{fit}{An object returned by a fitting function in this package.}
}
\details{
For standard joint Gauss-Hermite rules, the nodes are counted over the full integration dimension. For two-tier rules, each inner item rule is applied at each outer node. A subject may contribute fewer inner rules when an item is missing on all of that subject's occasions.
}
\value{
A list with the integration \code{method}, \code{layout}, \code{joint_dimension}, \code{largest_rule_dimension}, \code{largest_mode_dimension} and a data frame of \code{rules}. The joint dimension counts all random coordinates in the fitted representation. The largest mode dimension is zero for nonadaptive rules. Each row of \code{rules} gives the block name, dimension, number of copies and number of nodes per rule.
}
\seealso{\code{\link{growth_tt_fit}}, \code{\link{cfa_long_fit_md}}, \code{\link{cfa_long_fit_wu2016}}}
