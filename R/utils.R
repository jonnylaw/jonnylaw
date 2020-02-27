#' Title
#'
#' @param name
#'
#' @return
#' @export
#'
#' @examples
rename_clean <- function(name) {
  tolower(
    stringr::str_replace_all(string = name, pattern = "\\s+|/|-|>", replacement = "_")
  )
}
