#' Title
#'
#' @param samples
#' @param iters
#' @param parameters
#'
#' @return
#' @export
#'
#' @examples
iters_to_cube <- function(samples, iters, parameters) {
  samples %>%
    dplyr::sample_n(iters) %>%
    dplyr::select(tidyselect::contains(parameters)) %>%
    plyr::alply(1, diag)
}

#' Title
#'
#' Perform the Kalman filter using draws from the posterior distribution
#'
#' @param ws
#' @param vs
#' @param model
#'
#' @return
#' @export
#'
#' @examples
fitted_kalman_filter <- function(ys, ws, vs, model) {
  filtering <- purrr::map2(vs, ws, ~ kalmanFilter(ys = y, f = model$FF, g = model$GG, v = .x, w = .y, m0 = model$m0, c0 = model$C0))

  ys_filtered <- filtering %>%
    purrr::map("mt") %>%
    purrr::map(~ model$FF %*% .)

  ys_mean <- ys_filtered %>%
    purrr::reduce(`+`) / length(vs)

  lower <- abind::abind(ys_filtered, along = 3) %>%
    apply(c(1, 2), quantile, prob = 0.05)

  upper <- abind::abind(ys_filtered, along = 3) %>%
    apply(c(1, 2), quantile, prob = 0.95)

  list(mean = ys_mean, lower = lower, upper = upper)
}

#' Perform Kalman Smoother on a single site
#'
#' @param site
#' @param model
#'
#' @return
#' @export
#'
#' @examples
iterpolate_site <- function(site, model) {
  y <- site$value

  filtered <- dlm::dlmFilter(y = y, mod = model)
  smoothed <- dlm::dlmSmooth(filtered)

  smoothed_var <- dlm::dlmSvd2var(u = smoothed$U.S, d = smoothed$D.S)
  diagonal_var <- purrr::map(smoothed_var, diag) %>% unlist()

  tibble::tibble(
    date = site$date,
    smoothed_mean = smoothed$s[-1],
    upper = qnorm(0.975, mean = smoothed$s[-1], sd = sqrt(diagonal_var[-1])),
    lower = qnorm(0.025, mean = smoothed$s[-1], sd = sqrt(diagonal_var[-1])),
    y = y
  )
}
