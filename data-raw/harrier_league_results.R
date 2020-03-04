## code to prepare `harrier_league_results` dataset goes here
fixtures <- purrr::map_df(2013:2019, get_course_list)

# Get the raw html/xml document corresponding to each fixture
html_sen_men <- fixtures %>% 
  dplyr::rowwise() %>% 
  dplyr::filter(matching_course(course)) %>% 
  dplyr::mutate(course_formatted = format_course(course)) %>% 
  get_results_home(result_type = "SenM") %>% 
  dplyr::filter(any(class(results_home) == "xml_document")) %>% 
  dplyr::mutate(date = extract_date(results_home)) %>% 
  tibble::add_column(sex = "M")

html_sen_women <- fixtures %>% 
  dplyr::rowwise() %>% 
  dplyr::filter(matching_course(course)) %>% 
  dplyr::mutate(course_formatted = format_course(course)) %>% 
  get_results_home(result_type = "SenF") %>% 
  dplyr::filter(any(class(results_home) == "xml_document")) %>% 
  dplyr::mutate(date = extract_date(results_home)) %>% 
  tibble::add_column(sex = "F")

# parse the raw html/xml document
raw_results_table <- html_sen_women %>% 
  dplyr::bind_rows(html_sen_men) %>% 
  dplyr::group_by(date, course_formatted, sex) %>% 
  dplyr::do(results_table = try(parse_results_table(.$results_home[[1]])))

harrier_league_results <- raw_results_table %>% 
  tidyr::drop_na() %>% 
  dplyr::group_by(date, course_formatted, sex) %>%
  dplyr::do(results = clean_results(.$results_table[[1]])) %>% 
  tidyr::unnest() %>% 
  dplyr::mutate(
    time = ifelse(is.na(time), actual_time, time),
    time_seconds = as.numeric(lubridate::ms(time))
  )

usethis::use_data(harrier_league_results, overwrite = TRUE)
