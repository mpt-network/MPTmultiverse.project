#' Multiverse Analysis for MPT Models
#'
#' Performs a multiverse analysis for multinomial processing tree (MPT) models
#' across maximum-likelihood/frequentist and Bayesian estimation approaches. For
#' the frequentist approaches, no pooling (with and without parametric or
#' nonparametric bootstrap) and complete pooling  are implemented using
#' \pkg{MPTinR}. For the Bayesian approaches, no pooling, complete pooling, and
#' three different variants of partial pooling are implemented using
#' \pkg{TreeBUGS}. Requires \code{data} on a by-participant level with each row
#' corresponding to data from one participant (i.e., different response
#' categories correspond to different columns) and the data can contain a single
#' between-subjects condition. Model equations need to be passed as a
#' \code{.eqn} model file and category labels (first column in \code{.eqn} file)
#' need to match the column names in \code{data}. Results are returned in one
#' \code{tibble} with one row per estimation method.
#'
#' @param method \code{character} vector specifying which analysis approaches
#'   should be performed (see Description below). Defaults to all available
#'   methods.
#' @param dataset scalar \code{character} vector. Name of the data set that will
#'   be copied to the results \code{tibble}.
#' @param data A \code{data.frame} containing the data. Column
#'   names need to match category names in \code{model} (i.e., different from
#'   \pkg{MPTinR} behavior, order of categories is not important, matching is
#'   done via name).
#' @param model A model definition, typically the path to an \code{.eqn} model
#'   file containing the model equations. Category names need to match column
#'   names in \code{data}.
#' @param id scalar \code{character} vector. Name of the column that contains
#'   the subject identifier. If not specified, it is assumed that each row
#'   represents observations from one participant.
#' @param condition scalar \code{character} vector. Name of the column
#'   specifying a between-subjects factor. If not specified, no between-subjects
#'   comparisons are performed.
#' @param core \code{character} vector defining the core parameters of interest,
#'   e.g., \code{core = c("Dn", "Do")}. All other parameters are treated as
#'   auxiliary parameters.
#' @example examples/examples.fit_mpt.R
#' 
#' @details This functions is a fancy wrapper for packages \pkg{MPTinR} and
#'   \pkg{TreeBUGS} applying various frequentist and Bayesian estimation methods
#'   to the same data set using a single MPT model and collecting the results
#'   in one \code{tibble} where each row corresponds to one
#'   estimation method. Note that parameter restrictions (e.g., equating
#'   different parameters or fixing them to a constant) need to be part of the
#'   model (i.e., the \code{.eqn} file) and cannot be passed as an argument.
#'   
#'   The settings for the various methods are specified via function
#'   \code{\link{mpt_options}}. The default settings use all available cores for
#'   calculating the boostrap distribution as well as independent MCMC chains
#'   and should be appropriate for most situations.
#'
#'   The data can have a single between-subjects condition (specified via
#'   \code{condition}). This condition can have more than two levels. If
#'   specified, the pairwise differences between each level, the standard error
#'   of the differences, and confidence-intervals of the differences are
#'   calculated for each parameter. Please note that \code{condition} is
#'   silently converted to \code{character} in the output. Thus, a specific
#'   ordering of the \code{factor} levels in the output cannot be guaranteed.
#'   
#'   Parameter differences or other support for within-subject conditions is not
#'   provided. The best course of action for within-subjects conditions is to
#'   simply include separate trees and separate sets of parameters for each
#'   within-subjects condition. This allows to at least compare the estimates
#'   for each within-subjects condition across estimation method.
#'   
#'   \subsection{Implemented Methods}{
#'     Maximum-likelihood estimation with \pkg{MPTinR} via
#'     \code{\link[MPTinR]{fit.mpt}}:
#'     \itemize{
#'       \item{\code{"asymptotic_complete"}: }{Asymptotic ML theory, complete
#'       pooling}
#'       \item{\code{"asymptotic_no"}: }{ Asymptotic ML theory, no pooling}
#'       \item{\code{"pb_no"}: }{Parametric bootstrap, no pooling}
#'       \item{\code{"npb_no"}: }{Nonparametric bootstrap, no pooling}
#'     }
#'     
#'     Maximum-likelihood estimation with \pkg{HMMTreeR}
#'     \itemize{
#'       \item{\code{"latent_class"}: } { Asymptotic ML theory, partial pooling, latent-class approach}
#'     }
#'      
#'     Bayesian estimation with \pkg{TreeBUGS}
#'     \itemize{
#'       \item{\code{"simple"}: }{Bayesian estimation, no pooling (C++,
#'         \link[TreeBUGS]{simpleMPT})}
#'       \item{\code{"simple_pooling"}: }{Bayesian estimation, complete pooling 
#'         (C++, \link[TreeBUGS]{simpleMPT})}
#'       \item{\code{"trait"}: }{latent-trait model, partial pooling (JAGS, 
#'         \link[TreeBUGS]{traitMPT})}
#'       \item{\code{"trait_uncorrelated"}: }{latent-trait model without
#'         correlation parameters, partial pooling (JAGS,
#'         \link[TreeBUGS]{traitMPT})}
#'       \item{\code{"beta"}: }{beta-MPT model, partial pooling (JAGS, 
#'         \link[TreeBUGS]{betaMPT})}
#'       \item{\code{"betacpp"}: }{beta-MPT model, partial pooling (C++, 
#'         \link[TreeBUGS]{betaMPTcpp})}
#'     }
#'   }
#'   \subsection{Frequentist/Maximum-Likelihood Methods}{
#'     For the \emph{complete pooling asymptotic approach}, the group-level parameter
#'     estimates and goodness-of-fit statistics are the maximum-likelihood and
#'     G-squared values returned by \code{MPTinR}. The parameter differences are
#'     based on these values, the standard errors of the difference is simply
#'     the pooled standard error of the individual parameters. The overall fit
#'     (column \code{gof}) is based on an additional fit to the completely
#'     aggregated data.
#'     
#'     For the \emph{no pooling asymptotic approach}, the individual-level
#'     maximum-likelihood estimates are reported in column \code{est_indiv} and
#'     \code{gof_indiv} and provide the basis for the other results. Whether or
#'     not an individual-level parameter estimate is judged as identifiable
#'     (column \code{identifiable}) is based on separate fits with different
#'     random starting values. If, in these separate, fits the same objective
#'     criterion is reached several times (i.e., \code{Log.Likelihood} within
#'     .01 of best fit), but the parameter estimate differs (i.e., different
#'     estimates within .01 of each other), then an estimate is flagged as
#'     non-identifiable. If they are the same (i.e., within .01 of each other)
#'     they are marked as identifiable. The group-level parameters are simply
#'     the means of the identifiable individual-level parameters, the SE is the
#'     SE of the mean for these parameter (i.e., SD/sqrt(N), where N excludes
#'     non-identifiable parameters and thise estimated as NA), and the CI is
#'     based on mean and SE. The group-level and overall fit is the sum of the
#'     individual G-squares, sum of individual-level df, and corresponding
#'     chi-square df. The difference between the conditions and corresponding
#'     statistics are based on a t-test comparing the individual-level estimates
#'     (again, after excluding non-identifiable estimates). The CIs of the
#'     difference are based on the SEs (which are derived from a linear model
#'     equivalent to the t-test).
#'     
#'     The individual-level estimates of the \code{bootstrap based no-pooling}
#'     approaches are identical to the asymptotic ones. However, the SE is the
#'     SD of the bootstrapped distribution of parameter estimates, the CIs are
#'     the corresponding quantiles of the bootstrapped distribution, and the
#'     p-value is obtained from the bootstrapped G-square distribution.
#'     Identifiability of individual-level parameter estimates is also based on
#'     the bootstrap distribution of estimates. Specifically, we calculate the
#'     range of the CI (i.e., maximum minus minimum CI value) and flag those
#'     parameters as non-identifiable for which the range is larger than
#'     \code{mpt_options()$max_ci_indiv}, which defaults to \code{0.99}. Thus,
#'     in the default settings we say a parameter is non-identifiable if the
#'     bootstrap based CI extends from 0 to 1. The group-level estimates are the
#'     mean of the identifiable individual-level estimates. And difference
#'     between conditions is calculated in the same manner as for the asymptotic
#'     case using the identifiable individual-level parameter esatimates.
#'     
#'   }
#'   
#'   \subsection{Latent-class approach}{
#'     The \emph{latent-class approach} is fitted by interfacing \code{HMMTree},
#'     a software package that is only available on Microsoft Windows Machines.
#'     To install this software and the necessary \code{R} interface, use
#'     \code{devtools::install_github("methexp/HMMTreeR")}.
#'     It is currently not possible to estimate models that contain parameters
#'     that are fixed to numerical values.
#'     Multiple latent-class models with differing number of latent classes are
#'     estimated. The model that obtains the lowest \emph{AIC} while still being
#'     identified is selected for extracting parameter estimates
#'     The returned group-level parameter estimates are calculated as the
#'     weighted mean of parameter estimates of latent classes. Corresponding SEs
#'     are given by the square root of the weighted mean of class-wise squared
#'     SEs. Goodness-of-fit statistics are M1, M2, S1, and S2 as described by
#'     Klauer(2006).
#'   }
#'
#'   \subsection{Bayesian Methods}{
#'     The \emph{simple approaches} fit fixed-effects MPT models.
#'     \code{"simple"} uses no pooling and thus assumes independent uniform priors
#'     for the individual-level parameters. Group-level means are
#'     obtained as generated quantities by averaging the posterior samples
#'     across participants. \code{"simple_pooling"} aggregates observed
#'     frequencies across participants and assumes a uniform prior for the
#'     group-level parameters.
#'
#'     The \emph{latent-trait approaches} transform the individual-level
#'     parameters to a latent probit scale using the inverse cumulative standard
#'     normal distribution. For these probit values, a multivariate normal
#'     distribution is assumed at the group level. Whereas \code{"trait"}
#'     estimates the corresponding correlation matrix of the parameters
#'     (reported in the column \code{est_rho}), \code{"trait_uncorrelated"}
#'     assumes that the parameters are uncorrelated.
#'     
#'     For all Bayesian methods, the posterior distribution of the parameters is
#'     summarized by the posterior mean (in the column \code{est}), posterior
#'     standard deviation (\code{se}), and credbility intervals (\code{ci_*}).
#'     For parameter differences (\code{test_between}) and correlations
#'     (\code{est_rho}), Bayesian p-values are computed (column \code{p}) by
#'     counting the relative proportion of posterior samples that are smaller
#'     than zero. Goodness of fit is tested with the T1 statistic
#'     (observed vs. posterior-predicted average frequencies, \code{focus =
#'     "mean"}) and the T2 statistic (observed vs. posterior-predicted
#'     covariance of frequencies, \code{focus = "cov"}).
#'    }
#' 
#' @return A \code{tibble} with one row per estimation \code{method} and the
#'   following columns:
#' \enumerate{
#'   \item \code{model}: Name of model file (copied from \code{model} argument),
#'   \code{character}
#'   \item \code{dataset}: Name of data set (copied from \code{dataset}
#'   argument), \code{character}
#'   \item \code{pooling}: \code{character} specifying the level of pooling with
#'   three potential values: \code{c("complete", "no", "partial")}
#'   \item \code{package}: \code{character} specifying the package used for
#'   estimation with two potential values: \code{c("MPTinR", "TreeBUGS")}
#'   \item \code{method}: \code{character} specifying the method used with the
#'   following potential values: \code{c("asymptotic", "PB/MLE", "NPB/MLE",
#'   "simple", "trait", "trait_uncorrelated", "beta", "betacpp")}
#'   \item \code{est_group}: Group-level parameter estimates per condition/group.
#'   \item \code{est_indiv}: Individual-level parameter estimates (if provided
#'   by method). 
#'   \item \code{est_rho}: Estimated correlation of individual-level parameters
#'   on the probit scale (only in \code{method="trait"}). 
#'   \item \code{test_between}: Parameter differences between the levels of the
#'   between-subjects condition (if specified).
#'   \item \code{gof}: Overall goodness of fit across all individuals.
#'   \item \code{gof_group}: Group-level goodness of fit.
#'   \item \code{gof_indiv}: Individual-level goodness of fit.
#'   \item \code{fungibility}:  Posterior correlation of the group-level means 
#'   \code{pnorm(mu)} (only in \code{method="trait"}). 
#'   \item \code{test_homogeneity}: Chi-square based test of participant
#'   homogeneity proposed by Smith and Batchelder (2008). This test is the same
#'   for each estimation method.
#'   \item \code{convergence}: Convergence information provided by the
#'   respective estimation method. For the asymptotic frequentist methods this
#'   is a \code{tibble} with rank of the Fisher matrix, the number of parameters
#'   (which should match the rank of the Fisgher matrix), and the convergence
#'   code provided by the optimization algorithm (which is
#'   \code{\link{nlminb}}). The boostrap methods contain an additional column,
#'   \code{parameter}, that contains the information which (if any) parameters
#'   are empirically non-identifiable based on the bootstrapped distribution of
#'   parameter estimates (see above for exact description). For the Bayesian
#'   methods this is a \code{tibble} containing information of the posterior
#'   dsitribution (i.e., mean, quantiles, SD, SE, \code{n.eff}, and R-hat) for
#'   each parameter.
#'   \item \code{estimation}: Time it took for each estimation method and group.
#'   \item \code{options}: Options used for estimation. Obtained by running
#'   \code{\link{mpt_options}()}
#' }
#' 
#' With the exception of the first five columns (i.e., after \code{method}) all
#' columns are \code{list} columns typically holding one \code{tibble} per cell.
#' The simplest way to analyze the results is separately per column using
#' \code{\link[tidyr]{unnest}}. Examples for this are given below.
#'
#' @references 
#'   Smith, J. B., & Batchelder, W. H. (2008). Assessing individual differences
#'   in categorical data. \emph{Psychonomic Bulletin & Review}, 15(4), 713-731.
#'   \url{https://doi.org/10.3758/PBR.15.4.713}
#'   
#'   Klauer, K.C. (2006). Hierarchical multinomial processing tree models: A
#'   latent-class approach. \emph{Psychometrika}, \emph{71} (1), 7-31.
#'   \url{https://doi.org/10.1007/s11336-004-1188-3}
#'   
#' @importFrom MPTmultiverse fit_mpt
#' @export

fit_mpt <- function(
  model
  , dataset
  , data
  , id = NULL
  , condition = NULL
  , core = NULL
  , method
) {
  
  available_methods <- c(
    # MPTinR ----
    "asymptotic_complete"
    , "asymptotic_no"
    , "pb_no"
    , "npb_no"
    # TreeBUGS ----
    , "simple"
    , "simple_pooling"
    , "trait"
    , "trait_uncorrelated"
    , "beta"
    , "betacpp"
    # HMMTree ---
    , "latent_class"
  )
  
  if(missing(method)) {
    method <- available_methods
  }

  results <- MPTmultiverse::fit_mpt(
    model = model
    , dataset = dataset
    , data = data
    , id = id
    , condition = condition
    , core = core
    , method = union("asymptotic_complete", setdiff(method, "latent_class"))
  )

  results_lc <- list()
  # HMMTreeR part ----
  if("latent_class" %in% method) {

    # HMMTreeR installed?
    if(suppressWarnings(requireNamespace("HMMTreeR", quietly = TRUE))) {

      running_on_windows <- Sys.info()[["sysname"]]=="Windows"

      if(running_on_windows) {
        results_lc <- try(
          hmpt:::fit_lc(
            dataset = dataset
            , data = attr(results, "data")
            , model = attr(results, "model_file")
            , id = attr(results, "id")
            , condition = attr(results, "condition")
            , core = attr(results, "core")
          )
        )
      } else {
        message("Latent-class multinomial models can currently only be estimated on Windows -- sorry.")
      }
    } else {
      message("HMMTreeR is not installed on your system. Therefore, latent-class analyses are skipped.")
    }
  }

  if(tibble::is_tibble(results_lc)) {
    res <- list(
      multiverse = results
      , lc = results_lc
    )
    y <- dplyr::bind_rows(res)
    class(y) <- class(results)
    attr(y, "call") <- attr(results, "call")
    attr(y, "model_file") <- attr(results, "model_file")
    attr(y, "data_file") <- attr(results, "data_file")
    attr(y, "model") <- attr(results, "model")
    attr(y, "data") <- attr(results, "data")
    attr(y, "id") <- attr(results, "id")
    attr(y, "condition") <- attr(results, "condition")
    attr(y, "core") <- attr(results, "core")
  } else {
    y <- results
  }

  y
}