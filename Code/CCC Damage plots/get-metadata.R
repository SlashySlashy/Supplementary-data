library(tidyverse)
library(readxl)
get_metadata <- function(file) {
  cdata_file <- file
  med_cdata_basic <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_basic") %>%
    janitor::clean_names()
  med_cdata_extra <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_extra") %>%
    janitor::clean_names()
  med_cdata_cgg <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_cgg") %>%
    janitor::clean_names()
  med_cdata_dating <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_dating") %>%
    janitor::clean_names()
  med_cdata_geochem <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_geochemistry") %>%
    janitor::clean_names(replace = c("\u00b5" = "ug"))
  med_cdata_palinology <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_palinology") %>%
    janitor::clean_names()
  med_cdata_lipid_bulk <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_lipid_bulk") %>%
    janitor::clean_names()
  med_cdata_lipid_detailed <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_lipid_detailed") %>%
    janitor::clean_names()
  med_cdata_paleo_proxies <- readxl::read_xlsx(path = cdata_file, sheet = "med_metadata_paleoproxies") %>%
    janitor::clean_names()
  
  # We will keep the ones done in 2022, check cherry-picking.R for more details
  cherry_picked <- med_cdata_cgg %>%
    add_count(internal_name_by_dom, name = "n_reps") %>%
    filter(n_reps > 1, !is.na(sample_id)) %>%
    select(label = unique_sample_id, cgg_id) %>%
    filter(grepl("2021", label)) %>%
    pull(label)
  
  # med_cdata_cgg <- med_cdata_cgg %>%
  #   filter(!(unique_sample_id %in% cherry_picked))
  
  return(list(
    med_cdata_basic = med_cdata_basic,
    med_cdata_extra = med_cdata_extra,
    med_cdata_cgg = med_cdata_cgg,
    med_cdata_dating = med_cdata_dating,
    med_cdata_geochem = med_cdata_geochem,
    med_cdata_palinology = med_cdata_palinology,
    med_cdata_lipid_bulk = med_cdata_lipid_bulk,
    med_cdata_lipid_detailed = med_cdata_lipid_detailed,
    med_cdata_paleo_proxies = med_cdata_paleo_proxies
  ))
}
