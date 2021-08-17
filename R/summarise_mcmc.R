#' Thin MCMC
#'
#' Thin a dataframe by only keeping every nth row
#'
#' @param df a dataframe containing one sample per from from an MCMC
#' @param nth the number to thin by
#'
#' @return
#' @export
#'
#' @examples
thin <- function(df, nth) {
  df[seq(from = 1, to = nrow(df), by = nth), ]
}

#' Plot a traceplot
#'
#' @param chains
#' @param paginate boolean indicating whether to paginate the plots
#'
#' @return
#' @export
#'
#' @examples
traceplot <- function(chains) {
  if ("chain" %in% names(chains)) {
    chains %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = iteration, y = value, colour = as.factor(chain)), alpha = 0.5) +
      ggplot2::facet_wrap(~parameter, scales = "free_y", strip.position = "right") +
      ggplot2::theme(legend.position = "none")
  } else {
    chains %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = iteration, y = value), alpha = 0.5) +
      ggplot2::facet_wrap(~parameter, scales = "free_y", strip.position = "right") +
      ggplot2::theme(legend.position = "none")
  }
}

#' Plot a density plot from MCMC output
#'
#' @param chains
#'
#' @return
#' @export
#'
#' @examples
density_plot <- function(chains) {
  if ("chain" %in% names(chains)) {
    chains %>%
      ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = value, fill = as.factor(chain)), alpha = 0.5) +
      ggplot2::facet_wrap(~parameter, scales = "free", strip.position = "right") +
      ggplot2::theme(legend.position = "none")
  } else {
    chains %>%
      ggplot2::ggplot() +
      ggplot2::geom_density(ggplot2::aes(x = value), alpha = 0.5) +
      ggplot2::facet_wrap(~parameter, scales = "free", strip.position = "right") +
      ggplot2::theme(legend.position = "none")
  }
}

#' Title
#'
#' @param chains
#'
#' @return
#' @export
#'
#' @examples
plot_diagnostics <- function(chains) {
  p1 <- traceplot(chains)

  p2 <- density_plot(chains)

  gridExtra::grid.arrange(p1, p2)
}

#' Calculate Effective Sample Size
#'
#' @param values 
#'
#' @return
#' @export
#'
#' @examples
effective_size = function(values) {
  n = length(values)
  floor(n / (1 + 2 * sum(acf(x = values, plot = F)$acf)))
}

#' Plot the prior and posterior
#'
#' @param chain 
#' @param parameters 
#' @param prior 
#' @param limits 
#'
#' @return
#' @export
#'
#' @examples
prior_posterior = function(chain, parameters, prior, limits) {
  chain %>%
    dplyr::filter(parameter %in% parameters) %>%
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_line(aes(y = ..density.., colour = 'Posterior'), stat = 'density') +
    ggplot2::stat_function(fun = prior, ggplot2::aes(colour = 'Prior'), xlim = limits) +
    ggplot2::theme(legend.position = c(0.25, 0.9), title = ggplot2::element_blank()) +
    ggplot2::facet_wrap(~parameter, ncol = 1)
}

#' Title
#'
#' @param chains
#'
#' @return
#' @export
#'
#' @examples
summary_table <- function(chains) {
  chains %>%
    dplyr::group_by(Parameter) %>%
    dplyr::summarise(
      mean = mean(value), median = median(value),
      upper = quantile(value, probs = 0.95),
      lower = quantile(value, probs = 0.05),
      ESS = effective_size(value),
    )
}

#' Title
#'
#' @param chains
#'
#' @return
#' @export
#'
#' @examples
latex_summary_table <- function(chains) {
  summary_table(chains) %>%
    knitr::kable(digits = 2, booktabs = T, format = "latex")
}

#' Plot a probability density function
#'
#' @param pdf
#' @param mode
#' @param scale
#' @param range
#' @param title
#'
#' @return
#' @export
#'
#' @examples
#' plot_pdf()
plot_pdf <- function(pdf, mode, scale, range = c(-10, 10), title) {
  x <- seq(range[1], range[2], length.out = 1000)
  density <- pdf(x, mode, scale)
  ggplot2::qplot(x = x, y = density, geom = "line", xlim = range, main = title)
}

#' Title
#'
#' @param chains
#' @param actual_values
#'
#' @return
#' @export
#'
#' @examples
latex_table_sim <- function(chains, actual_values) {
  chains %>%
    dplyr::drop_na() %>%
    summary_table() %>%
    dplyr::inner_join(actual_values, by = "Parameter") %>%
    knitr::kable(digits = 2, booktabs = T, format = "latex")
}

#' Plot MCMC diagnostics with known values
#'
#' @param chains output from an MCMC sampler
#' @param actual_values dataframe containing the values used to simulate synthetic data
#'
#' @return
#' @export
#'
#' @examples
plot_diagnostics_sim <- function(chains, actual_values) {
  p1 <- chains %>%
    dplyr::inner_join(actual_values) %>%
    traceplot() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = actual_value), linetype = "dashed")

  p2 <- chains %>%
    dplyr::inner_join(actual_values) %>%
    density_plot() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = actual_value), linetype = "dashed")

  gridExtra::grid.arrange(p1, p2)
}
