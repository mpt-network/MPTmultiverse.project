#' Compare Estimates from Different Methods
#'
#' Returns all the information necessary to easily compare parameter estimates
#' from different analysis methods.
#'
#' @param x A results object from one or multiple multiverse analyses.
#' @param ... Further arguments that may be passed.
#'
#' @rdname compare_methods
#' @export

compare_methods <- function(x, ...) {
  UseMethod("compare_methods")
}

#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom graphics plot
#' @rdname compare_methods
#' @export

compare_methods.multiverseMPT <- function(
  x
  , ...
) {

  args <- list(...)

  if(is.null(args$which)) {
    args$which <- "est_group"
  }

  results <- tidyr::unnest(results, .data[[args$which]])

  results$inter <- interaction(
    results$method
    , results$pooling
    , results$package
    , drop = TRUE
    , sep = " "
  )

  results <- results %>%
    dplyr::group_by(.data$model, .data$dataset, .data$inter, .data$condition, .data$parameter) %>%
    dplyr::mutate(within = seq_along(.data$parameter)) %>%
    dplyr::ungroup()

  # Calculate CCC of all pairwise combinations
  pairs <- utils::combn(sort(levels(results$inter)), 2)
  all_pars_list <- vector("list", ncol(pairs))

  for (i in seq_len(ncol(pairs))) {
    tmp_dat <- results %>%
      dplyr::filter(.data$inter %in% pairs[, i]) %>%
      dplyr::select(.data$model, .data$dataset, .data$condition, .data$within, .data$parameter, .data$est, .data$inter) %>%
      dplyr::group_by(.data$model, .data$dataset) %>%
      dplyr::mutate(n_conditions = length(unique(.data$condition))) %>%
      dplyr::ungroup() %>%
      tidyr::spread(key = .data$inter, value = .data$est)
    colnames(tmp_dat)[(ncol(tmp_dat)-1):ncol(tmp_dat)] <- c("x", "y")
    tmp_dat$method_x <- pairs[1, i]
    tmp_dat$method_y <- pairs[2, i]
    all_pars_list[[i]] <- tmp_dat
  }

  all_pars <- dplyr::bind_rows(all_pars_list)
  all_pars
}

#' Pairs plot
#'
#' Create pairs plots of parameter estimates from multiple studies.
#'
#' @param x An object returned from \code{compare_methods}.
#' @param parameter Character. Indicating the name of the parameter to be plotted.
#' @param ... Further arguments that may be passed.
#'
#' @importFrom DescTools CCC
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export

pairs_plot <- function(x, parameter, ...) {

  if(!missing(parameter)) {
    x <- x[x$parameter == parameter, ]
  }
  
  format_ccc <- function(x) {
    ifelse(
      round(x, digits = 3) == 1
      , "1"
      , paste0(".", strsplit(format(x, digits = 3, nsmall = 3), split = "[.]")[[1]][2])
    )
  }

  plot_text <- x %>%
    dplyr::group_by(.data$method_x, .data$method_y) %>%
    dplyr::summarise(ccc = format_ccc(DescTools::CCC(.data$x, .data$y, na.rm = TRUE)$rho.c$est))
  
  x %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(method_x ~ method_y) +
    ggplot2::geom_text(
      data = plot_text
      , ggplot2::aes(
        x = 0.2
        , y = 0.9
        , label = plot_text$ccc
      )
      , parse = TRUE
      , inherit.aes = FALSE
      , size = 4
    ) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::xlab("") +
    ggplot2::ylab("")
}

#' Plot Concordance Correlation Coefficient
#'
#' Plot concordance correlation coefficients (CCCs) for parameters from multiple studies
#'
#' @param x An object returned from \code{compare_methods}.
#' @param ... Further arguments that may be passed.

#' @importFrom DescTools CCC
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export

ccc_plot <- function(x, ...) {
  # plot of CCC across parameters and conditions
  all_ccc <- x %>%
    dplyr::group_by(.data$method_x, .data$method_y, .data$parameter) %>%
    dplyr::summarise(
      ccc = DescTools::CCC(.data$x, .data$y, na.rm = TRUE)$rho.c$est
      , n_estimates = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    tibble::rowid_to_column("id")

  all_ccc %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$parameter, y = .data$ccc)) +
    ggbeeswarm::geom_beeswarm(ggplot2::aes_(shape = ~ method_x, col = ~ n_estimates),
                              cex = 0.5, alpha = 0.4) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se) +
    ggplot2::scale_colour_continuous(low = "#56B1F7", high = "#132B43") +
    ggplot2::scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13)) + 
    ggplot2::xlab("Parameter") + 
    ggplot2::ylab("CCC")
}
