#' Title
#'
#' @param name
#'
#' @return
#' @export
#'
#' @examples
rename_clean <- function(name) {
  janitor::make_clean_names(string = name)
}
