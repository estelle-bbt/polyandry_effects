#' Read data about plant flowers
#'
#' @description 
#' This function reads the data about plant flowers and format the table.
#'
#' @param file a character of length 1. The path to the .txt file.
#'
#' @return A `table` containing data. 
#' 
#' @import dplyr
#' 
#' @export

load_data_flower <- function(file_path)){
  
  data_flower <- read.table(file_path,head=T) |>
    mutate(poll_treat_factor=as.factor(case_when(poll_treat==1~"low",
                                                 poll_treat==2~"medium",
                                                 TRUE~"high"))) |>
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
  
  return(data_flower)
}

#' Read data about both sexual functions
#'
#' @description 
#'
#' @param 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

load_data_both_sexes <- function(file_path, c = "10")){
  
  data_both_sexes <- read.table(file_path,head=T) |>
    mutate(type=as.factor(type),
           poll_treat_factor=as.factor(poll_treat_factor)) |>
    mutate(type=relevel(type, ref = "fem"),
           poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high"))) |>
    mutate(loss_MS = 100-(100*gMS/get(paste0("nb_part_ID_out_co",c))))
  
  return(data_both_sexes)
}

#' Read data at the individual level
#'
#' @description 
#'
#' @param 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

load_data_id <- function(file_path)){
  
  data_id <- read.table(file_path,head=T) %>%
    mutate(poll_treat_factor=as.factor(case_when(poll_treat==1~"low",
                                                 poll_treat==2~"medium",
                                                 TRUE~"high"))) %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high"))) %>%
    mutate(prop_self=SR_self/(SR_self+SR_fem_out+SR_mal_out),
           gam_ov_proxy=nb_flo*nbOv_mean,
           SR_out=SR_fem_out+SR_mal_out) %>%
    mutate(mean_nb_visit_per_flower=nb_visit/nb_dist_vis)
  
  return(data_id)
}

#' Get data proxies of female reproductive success
#'
#' @description 
#'
#' @param 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

get_data_proxies <- function(data_flower, data_both_sexes, data_id, c = "10"){
  
  # resume information at the ID level to avoid pseudoreplication
  # we have to only retain row for which we have all informations according to each analyses
  # to avoid estimating each variable on a different number of sample for each individual
  
  # to check results without female that have less than 10 genotyped seeds
  # list_filt <- dti %>%
  #   filter(nGenot>=10) %>%
  #   pull(ID_full)
  
  data_proxy_ovav <- data_flower %>%
    # filter(ID_full %in% list_filt) %>%
    filter(!(is.na(nbOv)|is.na(nbAv)|is.na(nbGr_SR))) %>%
    group_by(ID_full,session,poll_treat_factor) %>%
    summarise(nbOv_sum=sum(nbOv,na.rm=T),
              nbAv_sum=sum(nbAv,na.rm=T),
              nbGr_sum=sum(nbGr_SR,na.rm=T)) %>%
    left_join(data_both_sexes %>% # data_both_sexes is the table with variables correctly implement as 0 or NA
                filter(type=="fem") %>%
                select(ID_full,gMS,paste0("nb_part_ID_out_co",c)) %>%
                rename(oMS=paste0("nb_part_ID_out_co",c))) %>%
    left_join(dti %>% select(ID_full, prop_self)) 
  
  data_proxy_ov <- data_flower %>%
    # filter(ID_full %in% list_filt) %>%
    filter(!(is.na(nbOv)|is.na(nbGr))) %>%
    group_by(ID_full,session,poll_treat_factor) %>%
    summarise(nbOv_sum=sum(nbOv,na.rm=T),
              nbGr_sum=sum(nbGr_SR,na.rm=T)) %>%
    left_join(data_both_sexes %>% # data_both_sexes is the table with variables correctly implement as 0 or NA
                filter(type=="fem") %>%
                select(ID_full,gMS,paste0("nb_part_ID_out_co",c)) %>%
                rename(oMS=paste0("nb_part_ID_out_co",c))) %>%
    left_join(dti %>% select(ID_full, prop_self))
  
  data_proxy_av <- data_flower %>%
    # filter(ID_full %in% list_filt) %>%
    filter(!(is.na(nbAv)|is.na(nbGr))) %>%
    group_by(ID_full,session,poll_treat_factor) %>%
    summarise(nbAv_sum=sum(nbAv,na.rm=T),
              nbGr_sum=sum(nbGr_SR,na.rm=T)) %>%
    left_join(data_both_sexes %>% # data_both_sexes is the table with variables correctly implement as 0 or NA
                filter(type=="fem") %>%
                select(ID_full,gMS,paste0("nb_part_ID_out_co",c)) %>%
                rename(oMS=paste0("nb_part_ID_out_co",c))) %>%
    left_join(dti %>% select(ID_full, prop_self))
  
  
  data_id_fem <- data_both_sexes %>%
    # filter(ID_full %in% list_filt) %>%
    filter(type=="fem") %>%
    rename(oMS=!!sym(paste0("nb_part_ID_out_co",c))) %>%
    left_join(data_id %>% select(ID_full, SR_fem_out,SR_self,seed_germ_sum,seed_semis_sum,poids_sd_mean,nb_frt,nb_flo,prop_self)) 
  
  return(list(data_proxy_ovav = data_proxy_ovav,
              data_proxy_ov = data_proxy_ov,
              data_proxy_av = data_proxy_av,
              data_id_fem = data_id_fem))
}