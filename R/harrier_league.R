#' Title
#'
#' @param raw_results
#'
#' @return
#'
#' @examples
clean_results <- function(raw_results) {
  raw_results %>%
    tibble::as_tibble() %>%
    dplyr::rename_all(janitor::make_clean_names) %>%
    dplyr::filter(cat != "guest") %>%
    dplyr::mutate(
      division = as.numeric(stringr::str_sub(club, start = 2, end = 2)),
      club = stringr::str_sub(club, start = 5)
    ) %>%
    dplyr::select(name, time, division, club) %>%
    dplyr::mutate(actual_seconds = as.numeric(lubridate::ms(time)))
}


#' Title
#'
#' @param url the url of a results page
#'
#' @return
#' @export
#'
#' @examples
get_harrier_league_results <- function(url) {
  htmltab::htmltab(doc = url) %>%
    clean_results()
}

#' Title
#'
#' @param strings
#'
#' @return
#'
#' @examples
concat <- function(strings) {
  purrr::reduce(.x = strings, .f = ~ paste(.x, .y, sep = ", "))
}

#' Get Harrier League Results
#'
#' This function calculates the position of a club in their respective division
#' by summing the positions of the first n = `counters` and sorting the sum
#' from smallest to largest.
#'
#' @param raw_results the finishing positions of the runners
#' @param counters the total number of counters to use for team results
#'
#' @return
#' @export
#'
#' @examples
get_division_results <- function(raw_results, counters = 4) {
  incomplete_teams = raw_results %>%
    count(club) %>%
    filter(n < counters)

  raw_results %>%
    dplyr::anti_join(incomplete_teams, by = "club") %>%
    dplyr::group_by(club) %>%
    dplyr::top_n(n = counters, wt = desc(pos)) %>%
    dplyr::group_by(club) %>%
    dplyr::summarise(
      total_points = sum(pos),
      total_counters = dplyr::n(),
      counters = concat(paste0(Name, " (", pos, ")"))
    ) %>%
    dplyr::mutate(position = min_rank(total_points))
}
