#' Get a list of courses for a given year
#'
#' @param year 
#'
#' @return
#' @export
#'
#' @examples
get_course_list <- function(year) {
  formatted_year <- paste(year, year + 1 - 2000, sep = "-")
  base_url <- "http://harrierleague.com/results/"
  year_url <- paste0(base_url, formatted_year)
  
  course_home <- xml2::read_html(year_url)
  
  result_urls <- course_home %>% 
    rvest::html_nodes(css = "body a") %>% 
    rvest::html_attr(name = "href")
  
  course_home %>% 
    rvest::html_nodes(css = "body a h3") %>% 
    rvest::html_text() %>% 
    tibble::tibble(course = .) %>% 
    tibble::add_column(year = year) %>% 
    tibble::add_column(result_url = result_urls)
}

#' Title
#'
#' @param course 
#'
#' @return
#' @export
#'
#' @examples
matching_course <- function(course) {
  possible_courses <-
    c(
      "alnwick",
      "aykley heads",
      "cramlington",
      "druridge",
      "herrington",
      "jarrow",
      "gosforth",
      "sherman",
      "thornley",
      "wrekenton",
      "wallington"
    )
  any(purrr::map_lgl(possible_courses, grepl, x = course, ignore.case = T))
}

#' Title
#'
#' @param course 
#'
#' @return
#' @export
#'
#' @examples
format_course <- function(course) {
  possible_courses <-
    c(
      "alnwick",
      "aykley heads",
      "cramlington",
      "druridge",
      "herrington",
      "jarrow",
      "gosforth",
      "sherman",
      "thornley",
      "wrekenton",
      "wallington"
    )
  bools <- purrr::map_lgl(possible_courses, grepl, x = course, ignore.case = T)
  possible_courses[bools]
}

#' Title
#'
#' @param parsed_courses 
#' @param possible_courses 
#'
#' @return
#' @export
#'
#' @examples
format_course_name <- function(parsed_courses, possible_courses) {
  parsed_courses %>% 
    dplyr::rowwise() %>% 
    dplyr::filter(matching_course(course)) %>% 
    dplyr::mutate(course_formatted = format_course(course))
}

#' Get the URL of a single fixture
#'
#' @param courses 
#' @param result_type 
#'
#' @return
#' @export
#'
#' @examples
get_results_home <- function(courses, result_type = "SenM") {
  base_url <- "http://harrierleague.com"
  courses %>% 
    dplyr::group_by(year, course_formatted, result_url) %>% 
    dplyr::do(results_home = try(xml2::read_html(URLencode(paste0(base_url, .$result_url, result_type, ".htm"))), silent = TRUE))
}

#' Get the date of a fixture
#' 
#' Extract the date from an h1 or h2 title on the results page
#'
#' @param results_home the url of a fixture
#'
#' @return
#' @export
#'
#' @examples
extract_date <- function(results_home) {
  date <- results_home %>%
    rvest::html_nodes("body h1") %>%
    rvest::html_text() %>%
    concat() %>%
    stringr::str_extract(., "[0-9]{1,2}.([0-9]{2}|[A-Za-z]+).[0-9]{4}") %>%
    lubridate::parse_date_time(x = ., orders = c("d-m-Y", "m-d-Y", "d b Y"))
  
  if (is.na(date)) {
    date <- results_home %>%
      rvest::html_nodes("body h2") %>%
      rvest::html_text() %>%
      concat() %>%
      stringr::str_extract(., "[0-9]{1,2}.([0-9]{2}|[A-Za-z]+).[0-9]{4}") %>%
      lubridate::parse_date_time(x = ., orders = c("d-m-Y", "m-d-Y", "d b Y"))
  }
  
  date
}

#' Get the raw HTML table from HTML
#'
#' @param html a string containing html for the page
#'
#' @return
#' @export
#'
#' @examples
parse_results_table <- function(html) {
  html %>% 
    rvest::html_table() %>% 
    purrr::pluck(1)
}

#' Check if some columns contain non-missing values
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
not_all_na <- function(x) {
  any(!is.na(x))
}

#' Clean one fixture of scraped results
#'
#' @param raw_results
#'
#' @return
#'
#' @export
#' @examples
clean_results <- function(raw_results) {
  raw_results %>%
    dplyr::select_if(.predicate = not_all_na) %>% 
    tibble::as_tibble() %>%
    dplyr::rename_all(janitor::make_clean_names) %>%
    dplyr::filter(cat != "guest") %>%
    dplyr::mutate(
      name = stringr::str_extract(name, pattern = "[A-Za-z&\\s+]+"),
      division = stringr::str_extract(club, "[1-3]"),
      club = stringr::str_extract(club, pattern = "[A-Za-z&\\s+]+")
    ) %>%
    dplyr::select(name, cat, dplyr::contains("time"), division, club)
}


#' Fetch a table containing the harrier league results
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

#' Concatenate strings together as a comma-separated list
#'
#' @param strings a list of strings
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
  incomplete_teams <- raw_results %>%
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
