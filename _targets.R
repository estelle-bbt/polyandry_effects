#' Targets plan
#' 

## Attach required packages ----

library(targets)
library(tarchetypes)
library(ggplot2)

tar_option_set(
  packages = c("dplyr","tidyr","ggplot2","forcats")  # load dplyr in each environement
)

tar_source()

## Load Project R Functions ----

source(here::here("R", "functions.R"))

## Analyses pipeline ----

list(
  
  ## Manage data ----
  
  tar_target(data_flower,load_data_flower("data/data_ABPOLL_flower_resume.txt")),
  
  tar_target(data_id,load_data_id("data/data_ABPOLL_ID_resume.txt")),
  
  tar_target(data_both_sexes,load_data_both_sexes("data/all_data_long_NA_0AllFemFALSE_raw.txt")),
  
  tar_target(data_proxies,get_data_proxies(data_flower, data_both_sexes, data_id, c = "10")),
  
  tar_target(common_theme,get_common_theme()),
  
  tar_target(prop_non_abort,get_prop_non_abort(data_proxies$data_proxy_av, common_theme)),
  
  tar_target(seed_set,get_seed_set(data_proxies$data_proxy_ov, common_theme)),
  
  tar_target(prop_germ,get_prop_germ(data_proxies$data_id_fem, common_theme)),
  
  tar_target(weight,get_weight(data_proxies$data_id_fem, common_theme)),
  
  
  
  ## Quarto ----
  
  tarchetypes::tar_quarto(index, "index.qmd")
  
)



