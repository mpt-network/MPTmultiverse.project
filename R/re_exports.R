#' Options Settings for MPT Comparison
#' 
#' Set and examine a variety of \emph{options} which affect the way MPT models
#' are estimated.
#' 
#' @param ... Named parameters to set. Possible values are:
#' \itemize{
#'   \item{\code{bootstrap_samples}: }{Numeric. The number of bootstrap samples to be drawn for the calculation parametric bootstrap confidence intervals.}
#'   \item{\code{n.optim}: }{Numeric. The number of optimization runs for the models estimated with maximum-likelihood methods.}
#'   \item{\code{n.chains}: }{Numeric. The number of MCMC chains for the Bayesian models.}
#'   \item{\code{n.adapt}: }{Numeric. The number of iterations for adaptation.}
#'   \item{\code{n.burnin}: }{Numeric. The number of burn-in/warm-up iterations.}
#'   \item{\code{n.iter}: }{Numeric. The total number of iterations to be drawn \emph{after} adaptation (including burnin).}
#'   \item{\code{n.thin}: }{Numeric. Thinning interval.}
#'   \item{\code{Rhat_max}: }{Numeric. The maximum rhat.}
#'   \item{\code{Neff_min}: }{Numeric. The minimum number of effective samples you are willing to accept.}
#'   \item{\code{extend_max}: }{Numeric.}
#'   \item{\code{n.PPP}: }{Numeric. The number of posterior predictive samples drawn for the calculation of fit statistics T_1 and T_2.}
#'   \item{\code{n.CPU}: }{Numeric. The number of CPU cores to use for obtaining the parametric bootstrap dsitribution. Defaults to the number of available cores on your machine.}
#'   \item{\code{ci_size}: }{Numeric.}
#'   \item{\code{max_ci_indiv}: }{Numeric. Used for excluding individual parameter estimates in the bootstrap approaches. If the range of the CI (i.e., distance between minimum and maximum) is larger than this value, the estimate is excluded from the group-level estimates.}
#'   \item{\code{silent_jags}: }{Logical. Whether to suppress JAGS output.}
# ' TODO  \item{\code{catch_warnings}: }{Logical. Whether to store warnings and errors as additional columns in the output.}
#'   \item{\code{save_models}: }{Logical.}
#' }
#'   
#' 
#' @examples 
#' # Examine options:
#' mpt_options()
#' 
#' # Set number of MCMC chains to 20:
#' mpt_options(n.chains = 20)
#' mpt_options()
#' 
#' @importFrom MPTmultiverse mpt_options
#' @export
mpt_options <- MPTmultiverse::mpt_options


#' Check results from a multiverse analysis
#' 
#' This is a helper function to see if the model estimation worked as intended.
#' 
#' @param results An object of class multiverseMPT.
#' @importFrom MPTmultiverse check_results
#' @export
check_results <- MPTmultiverse::check_results


#' Write check_results
#'
#' Helper function to write the output from check_results() to a file.
#'
#' @param DATA_FILE Character. File name to use.
#' @param results An object of class multiverseMPT.
#' @importFrom MPTmultiverse write_check_results
#' @export
write_check_results <- MPTmultiverse::write_check_results
