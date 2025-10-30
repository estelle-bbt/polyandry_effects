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
  
  tar_target(file_flower,"data/data_ABPOLL_flower_resume.txt",format = "file"), # to reevaluate the next target if the file is modified
  
  tar_target(data_flower,load_data_flower(file_flower)),
  
  tar_target(oms_gms_id_flower,get_oms_gms_id_flower(file_path_obs = "data/data_ABPOLL_ID_level_detflo.txt", file_path_gen = "data/fix10_paternities_ABPOLL.txt")),
  
  tar_target(data_proxy_id, get_data_proxy_id(data_flower)),
  
  tar_target(data_final,get_data_final(oms_gms_id_flower$oms_gms_flower, oms_gms_id_flower$oms_gms_id, data_flower, data_proxy_id)),
  
  tar_target(stat_flower_wo_contact,get_stat_flower(data_final$data_final_flower, include_contact = FALSE)),

  tar_target(stat_flower_w_contact,get_stat_flower(data_final$data_final_flower, include_contact = TRUE)),
  
  tar_target(stat_id_wo_contact,get_stat_id(data_final$data_final_id, include_contact = FALSE)),
  
  tar_target(stat_id_w_contact,get_stat_id(data_final$data_final_id, include_contact = TRUE)),
  
  tar_target(results_plot,get_results_plot(stat_flower_wo_contact,stat_flower_w_contact,stat_id_wo_contact,stat_id_w_contact)),

  # tar_target(data_id,load_data_id("data/data_ABPOLL_ID_resume.txt")),
  
  # tar_target(data_both_sexes,load_data_both_sexes("data/all_data_long_NA_0AllFemFALSE_raw.txt")),
  
  # tar_target(data_proxies,get_data_proxies(data_flower, data_both_sexes, data_id, c = "10")),
  # 
  # tar_target(common_theme,get_common_theme()),
  # 
  # tar_target(prop_non_abort,get_prop_non_abort(data_proxies$data_proxy_av, common_theme)),
  # 
  # tar_target(seed_set,get_seed_set(data_proxies$data_proxy_ov, common_theme)),
  # 
  # tar_target(prop_germ,get_prop_germ(data_proxies$data_id_fem, common_theme)),
  # 
  # tar_target(weight,get_weight(data_proxies$data_id_fem, common_theme)),
  # 
  # 
  
  ## Quarto ----
  
  tarchetypes::tar_quarto(index, "index.qmd")
  
)



