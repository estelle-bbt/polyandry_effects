#' Targets plan
#' 

## Attach required packages ----

library(targets)
library(tarchetypes)
library(ggplot2)

tar_option_set(
  packages = c("dplyr","tidyr","ggplot2")  # load dplyr in each environement
)

tar_source()

## Load Project R Functions ----

source(here::here("R", "functions.R"))

## Analyses pipeline ----

list(
  
  ## Manage data ----
  
  tar_target(data_flower,load_data_flower("data/data_ABPOLL_flower_resume.txt")),
  
  tar_target(data_id,load_data_id("data/data_ABPOLL_ID_resume.txt")),
  
  tar_target(data_both_sexes,load_data_both_sexes("data/all_data_long_NA_0AllFemFALSE.txt")),
  
  tar_target(data_proxies,get_data_proxies(data_flower, data_both_sexes, data_id, c = "10")),
  
  
  ## Quarto ----
  
  tarchetypes::tar_quarto(index, "index.qmd")
  
)



