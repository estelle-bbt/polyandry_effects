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

load_data_flower <- function(file_path){
  
  data_flower <- read.table(file_path,head=T) |>
    filter(frt == 1) |> # we only keep the fruit for rs proxies
    mutate(ttt=as.factor(case_when(poll_treat==1~"low",
                                                 poll_treat==2~"medium",
                                                 TRUE~"high"))) |>
    mutate(ttt=forcats::fct_relevel(ttt, c("low","medium","high"))) |>
    rename(id = ID_full,
           nb_ab = nbAv,
           nb_ov = nbOv,
           nb_seed = nbGr_SR,
           seed_weight = poids_sd,
           nb_germ = nbGr_germ,
           nb_sown = nbGr_semis)
  
  return(data_flower)
}

#' #' Read data about both sexual functions
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export
#' 
#' load_data_both_sexes <- function(file_path, c = "10"){
#'   
#'   data_both_sexes <- read.table(file_path,head=T) |>
#'     mutate(type=as.factor(type),
#'            poll_treat_factor=as.factor(poll_treat_factor)) |>
#'     mutate(type=relevel(type, ref = "fem"),
#'            poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high"))) |>
#'     mutate(loss_MS = 100-(100*gMS/get(paste0("nb_part_ID_out_co",c))))
#'   
#'   return(data_both_sexes)
#' }

#' Estimate q index
#'
#' @description 
#'
#' @param data 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

q_index <- function(x) {
  n <- length(x)
  N <- sum(x)
  if (N <= 1) return(NA)  # to avoid dividing by zero
  sum_xi <- sum(x * (x - 1))
  I_d <- (sum_xi) / (N * (N - 1))
  return(I_d)
}

#' Get table with oms/gms at the id and flower scale
#'
#' @description filtering seeds that cannot be there
#'
#' @param 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

get_oms_gms_id_flower <- function(file_path_obs = "data/data_ABPOLL_ID_level_detflo.txt", file_path_gen = "data/fix10_paternities_ABPOLL.txt"){
  
  # flower that only received self pollen = true 0 MS
  only_self_flower <- read.table(file_path_obs,head=T) |>
    filter(!session %in% c("5.FA1","5.MO1"))  |>
  group_by(ID_full_part, id_flow_part) |>
  summarise(
    # nb individual with export
    n_nonzero = sum(export_nb_visit_co10 != 0),
    # si la seule ligne non nulle correspond à un "self"
    self_only = ifelse(
      n_nonzero == 1 &
      all(ID_full_foc[export_nb_visit_co10 != 0] == ID_full_part[export_nb_visit_co10 != 0]),
      1, 0
    ),
    .groups = "drop"
  ) |>
  select(ID_full_part, id_flow_part, self_only) |>
    filter (self_only == 1) |>
    mutate(oms = 0) |> 
    rename(id = ID_full_part,
           id_flow = id_flow_part) |>
    select(-self_only)
  
  # id that only received self pollen = true 0 MS
  only_self_id <- read.table(file_path_obs,head=T) |>
    filter(!session %in% c("5.FA1","5.MO1"))  |>
    group_by(ID_full_part) |>
    summarise(
      # nb individual with export
      n_nonzero = sum(export_nb_visit_co10 != 0),
      # si la seule ligne non nulle correspond à un "self"
      self_only = ifelse(
        n_nonzero == 1 &
          all(ID_full_foc[export_nb_visit_co10 != 0] == ID_full_part[export_nb_visit_co10 != 0]),
        1, 0
      ),
      .groups = "drop"
    ) |>
    select(ID_full_part, self_only) |>
    filter (self_only == 1) |>
    mutate(oms = 0) |> 
    rename(id = ID_full_part) |>
    select(-self_only)
  # not any id received only self pollen
    
  obs_clean <- read.table(file_path_obs,head=T) |>
    filter(!session %in% c("5.FA1","5.MO1"))  |>
    select(session,ID_full_foc,ID_full_part,id_flow_part,export_nb_visit_co10) |>
    filter(ID_full_foc != ID_full_part)  |>
    filter(!(is.na(export_nb_visit_co10)|export_nb_visit_co10==0)) |>
    mutate(couple_obs = paste0(ID_full_foc,"_",id_flow_part)) |>
    rename(id = ID_full_part,
           id_flow = id_flow_part)
  
  # filtering only observed contact
  gen_clean <- read.table(file_path_gen,head=T) |>
    filter(known_id != candidate_id) |>
    mutate(offspring_id = stringr::str_remove(offspring_id, pattern = "_run3"),
           offspring_id = stringr::str_remove(offspring_id, pattern = "a"),
           offspring_id = stringr::str_remove(offspring_id, pattern = "b")) |>
    mutate(flower_id = sub("_[^_]+$", "", offspring_id)) |>
    mutate(couple_gen = paste0(candidate_id,"_",flower_id)) |>
    filter(couple_gen %in% obs_clean$couple_obs)  |>
    rename(id = known_id,
           id_flow = flower_id)
  
  q_obs_id <- obs_clean |>
    group_by(id) %>%
    summarise(
      q_obs = q_index(export_nb_visit_co10),
      .groups = "drop"
    ) 
  
  q_gen_id <- gen_clean |>
    group_by(id,candidate_id) |>
    summarise(genot_couple = n()) |>
    group_by(id) |>
    summarise(
      q_gen = q_index(genot_couple),
      .groups = "drop"
    )
  
  q_obs_flower <- obs_clean |>
    group_by(id, id_flow) %>%
    summarise(
      q_obs = q_index(export_nb_visit_co10),
      .groups = "drop"
    ) 
  
  q_gen_flower <- gen_clean |>
    group_by(id, id_flow, candidate_id) |>
    summarise(genot_couple = n()) |>
    group_by(id, id_flow) |>
    summarise(
      q_gen = q_index(genot_couple),
      .groups = "drop"
    )
  
  oms_gms_flower <- obs_clean |> 
    group_by(id,id_flow) |>
    summarise(oms = n_distinct(ID_full_foc),
              contact = mean(export_nb_visit_co10, na.rm = T))  |>
    bind_rows(only_self_flower) |>
    left_join(q_obs_flower) |>
    left_join(q_gen_flower) |>
    mutate(ratio_q = q_gen / q_obs)
  
  oms_gms_id <- obs_clean |> 
    group_by(id) |>
    summarise(oms = n_distinct(ID_full_foc),
              contact = mean(export_nb_visit_co10, na.rm = T))  |>
    # bind_rows(only_self_id) |> # not any id
    left_join(q_obs_id) |>
    left_join(q_gen_id) |>
    mutate(ratio_q = q_gen / q_obs)
  
  return(list(oms_gms_flower = oms_gms_flower,
         oms_gms_id = oms_gms_id))
}

#' #' Get data proxies female rs at the individual level
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export

get_data_proxy_id <- function(data_flower){

  # resume information at the ID level to avoid pseudoreplication
  # we have to only retain row for which we have all informations according to each analyses
  # to avoid estimating each variable on a different number of sample for each individual

  # proportion non-aborted
  data_proxy_ab <- data_flower %>%
    filter(!(is.na(nb_ab)|is.na(nb_seed))) %>%
    group_by(id,session,ttt) %>%
    summarise(nb_ab_sum_ab=sum(nb_ab,na.rm=T),
              nb_seed_sum_ab=sum(nb_seed,na.rm=T))
  
  # seed-set
  data_proxy_ss <- data_flower %>%
    filter(!(is.na(nb_ov)|is.na(nb_seed))) %>%
    group_by(id,session,ttt) %>%
    summarise(nb_ov_sum_ss=sum(nb_ov,na.rm=T),
              nb_seed_sum_ss=sum(nb_seed,na.rm=T))
  
  # proportion germinated
  data_proxy_germ <- data_flower %>%
    filter(!(is.na(nb_sown)|is.na(nb_germ))) %>%
    group_by(id,session,ttt) %>%
    summarise(nb_sown_sum_germ=sum(nb_sown,na.rm=T),
              nb_germ_sum_germ=sum(nb_germ,na.rm=T))
  
  # proportion germinated
  data_proxy_weight <- data_flower %>%
    group_by(id,session,ttt) %>%
    summarise(mean_seed_weight=mean(seed_weight,na.rm=T))
  
  data_proxy_id <- data_proxy_ab |>
    full_join(data_proxy_ss) |>
    full_join(data_proxy_germ) |>
    full_join(data_proxy_weight) |>
    ungroup() |>
    mutate(across(4:10, ~ ifelse(is.nan(.), NA, .))) |> # for id without seed weight
    filter(if_any(4:10, ~ !is.na(.)))  # remove id without data
    
  return(data_proxy_id)
}

#' #' Get final data table with rs proxy and oms/gms
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export

get_data_final <- function(oms_gms_flower, oms_gms_id, data_flower, data_proxy_id){
  
  data_final_flower <- data_flower |> # only receptive flowers
    select(id, id_flow, session, ttt, nb_ab, nb_ov, nb_seed, seed_weight, nb_sown, nb_germ, pl) |>
    left_join(oms_gms_flower) |>
    filter(!is.na(oms)) # filter spontaneous selfing or obs error
  
  data_final_id <- data_proxy_id |>
    left_join(oms_gms_id) |>
    filter(!is.na(oms)) # filter spontaneous selfing or obs error
  
  return(list(data_final_flower = data_final_flower,
              data_final_id = data_final_id))
}

#' Variance intra vs inter pollen load
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

get_var_intra_inter <- function(data_final_flower){
  
  # intra-class correlation coefficient
  
  model_pl <- lme4::glmer(data = data_final_flower, pl ~ 1 + (1 | id) + (1|session), family = "poisson")
  
  model_pl <- lme4::glmer(data = data_final_flower, pl ~ 1 + (1|session) + (1 |session:id), family = "poisson")
  
  
  vc <- as.data.frame(lme4::VarCorr(model_pl))
  var_id <- vc[vc$grp == "session:id", "vcov"]
  var_session <- vc[vc$grp == "session", "vcov"]
  
  icc_id <- var_id / (var_id + var_session)
  # 0.77 not very far from 1 -> var intra < var inter
  # ok to resume pl with its mean at the individual level
  
  return(icc_id)
}

#' Statistic models: flower scale
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

get_stat_flower <- function(data_final_flower){
  
  # effect of oms on q ratio ----
  model_q_oms_flower_0 <- lme4::lmer(data = data_final_flower, ratio_q ~ oms * ttt + pl * ttt + (1|session) + (1|session:id))
  model_q_oms_flower_1 <- lme4::lmer(data = data_final_flower, ratio_q ~ oms + ttt + pl * ttt + (1|session) + (1|session:id))
  model_q_oms_flower_2 <- lmer4::lmer(data = data_final_flower, ratio_q ~ oms * ttt + pl + ttt + (1|session) + (1|session:id))
  model_q_oms_flower_3 <- lme4::lmer(data = data_final_flower, ratio_q ~ oms + ttt + pl + (1|session) + (1|session:id))
  model_q_oms_flower_3b <- lme4::lmer(data = data_final_flower   |> filter (!is.na(pl)), ratio_q ~ oms + ttt + pl + (1|session) + (1|session:id))
  model_q_oms_flower_4 <- lme4::lmer(data = data_final_flower, ratio_q ~ ttt + pl + (1|session) + (1|session:id))
  model_q_oms_flower_5 <- lme4::lmer(data = data_final_flower, ratio_q ~ oms + pl + (1|session) + (1|session:id))
  model_q_oms_flower_6 <- lme4::lmer(data = data_final_flower  |> filter (!is.na(pl)), ratio_q ~ oms + ttt + (1|session) + (1|session:id))
  anova(model_q_oms_flower_0, model_q_oms_flower_1) # sign oms * ttt
  anova(model_q_oms_flower_0, model_q_oms_flower_2) # ns pl * ttt
  anova(model_q_oms_flower_3, model_q_oms_flower_4) # sign oms
  anova(model_q_oms_flower_3, model_q_oms_flower_5) # marg sign ttt
  anova(model_q_oms_flower_3b, model_q_oms_flower_6) # ns pl
  
  data_final_flower <- data_final_flower |>
    mutate(ttt=forcats::fct_relevel(ttt, "high"))
  model_q_oms_flower_best <- lmerTest::lmer(data = data_final_flower, ratio_q ~ oms * ttt + pl + ttt + (1|session) + (1|session:id))
  summary(model_q_oms_flower_best)
  
  # effect of oms on rs proxies ----
  ## aborted ----
  model_ab_oms_flower_0 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ oms * ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ab_oms_flower_1 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ oms + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ab_oms_flower_2 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ oms * ttt + pl + ttt + (1|session) + (1|session:id), family = "binomial")
  model_ab_oms_flower_3 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ oms + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_oms_flower_3b <- lme4::glmer(data = data_final_flower  |> filter (!is.na(pl)), cbind(nb_seed,nb_ab) ~ oms + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_oms_flower_4 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_oms_flower_5 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ oms + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_oms_flower_6 <- lme4::glmer(data = data_final_flower  |> filter (!is.na(pl)), cbind(nb_seed,nb_ab) ~ oms + ttt + (1|session) + (1|session:id), family = "binomial")
  anova(model_ab_oms_flower_0, model_ab_oms_flower_1) # ns oms * ttt
  anova(model_ab_oms_flower_0, model_ab_oms_flower_2) # sign pl * ttt
  anova(model_ab_oms_flower_3, model_ab_oms_flower_4) # ns oms
  anova(model_ab_oms_flower_3, model_ab_oms_flower_5) # ns sign ttt
  anova(model_ab_oms_flower_3b, model_ab_oms_flower_6) # ns pl

  data_final_flower <- data_final_flower |>
    mutate(ttt=forcats::fct_relevel(ttt, "high"))
  model_ab_oms_flower_best <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ oms + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  summary(model_ab_oms_flower_best)
  
  ## seed-set ----
  model_ss_oms_flower_0 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ oms * ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ss_oms_flower_1 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ oms + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ss_oms_flower_2 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ oms * ttt + pl + ttt + (1|session) + (1|session:id), family = "binomial")
  model_ss_oms_flower_3 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ oms + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_oms_flower_3b <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_seed, nb_ov - nb_seed) ~ oms + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_oms_flower_4 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_oms_flower_5 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ oms + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_oms_flower_6 <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_seed, nb_ov - nb_seed) ~ oms + ttt + (1|session) + (1|session:id), family = "binomial")
  anova(model_ss_oms_flower_0, model_ss_oms_flower_1) # sign oms * ttt
  anova(model_ss_oms_flower_0, model_ss_oms_flower_2) # sign pl * ttt
  anova(model_ss_oms_flower_3, model_ss_oms_flower_4) # sign oms
  anova(model_ss_oms_flower_3, model_ss_oms_flower_5) # ns ttt
  anova(model_ss_oms_flower_3b, model_ss_oms_flower_6) # sign pl

  data_final_flower <- data_final_flower |>
    mutate(ttt=forcats::fct_relevel(ttt, "high"))
  model_ss_oms_flower_best <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ oms * ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  summary(model_ss_oms_flower_best)
  
  ## germ ----
  model_germ_oms_flower_0 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ oms * ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_germ_oms_flower_1 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ oms + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_germ_oms_flower_2 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ oms * ttt + pl + ttt + (1|session) + (1|session:id), family = "binomial")
  model_germ_oms_flower_3 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ oms + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_oms_flower_3b <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_germ, nb_sown - nb_germ) ~ oms + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_oms_flower_4 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_oms_flower_5 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ oms + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_oms_flower_6 <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_germ, nb_sown - nb_germ) ~ oms + ttt + (1|session) + (1|session:id), family = "binomial")
  anova(model_germ_oms_flower_0, model_germ_oms_flower_1) # ns oms * ttt
  anova(model_germ_oms_flower_0, model_germ_oms_flower_2) # ns pl * ttt
  anova(model_germ_oms_flower_3, model_germ_oms_flower_4) # ns oms
  anova(model_germ_oms_flower_3, model_germ_oms_flower_5) # ns ttt
  anova(model_germ_oms_flower_3b, model_germ_oms_flower_6) # ns pl

  ## weight ----
  model_weight_oms_flower_0 <- lme4::lmer(data = data_final_flower, seed_weight ~ oms * ttt + pl * ttt + (1|session) + (1|session:id))
  model_weight_oms_flower_1 <- lme4::lmer(data = data_final_flower, seed_weight ~ oms + ttt + pl * ttt + (1|session) + (1|session:id))
  model_weight_oms_flower_2 <- lme4::lmer(data = data_final_flower, seed_weight ~ oms * ttt + pl + ttt + (1|session) + (1|session:id))
  model_weight_oms_flower_3 <- lme4::lmer(data = data_final_flower, seed_weight ~ oms + ttt + pl + (1|session) + (1|session:id))
  model_weight_oms_flower_3b <- lme4::lmer(data = data_final_flower |> filter (!is.na(pl)), seed_weight ~ oms + ttt + pl + (1|session) + (1|session:id))
  model_weight_oms_flower_4 <- lme4::lmer(data = data_final_flower, seed_weight ~ ttt + pl + (1|session) + (1|session:id))
  model_weight_oms_flower_5 <- lme4::lmer(data = data_final_flower, seed_weight ~ oms + pl + (1|session) + (1|session:id))
  model_weight_oms_flower_6 <- lme4::lmer(data = data_final_flower |> filter (!is.na(pl)), seed_weight ~ oms + ttt + (1|session) + (1|session:id))
  anova(model_weight_oms_flower_0, model_weight_oms_flower_1) # ns oms * ttt
  anova(model_weight_oms_flower_0, model_weight_oms_flower_2) # ns pl * ttt
  anova(model_weight_oms_flower_3, model_weight_oms_flower_4) # ns oms
  anova(model_weight_oms_flower_3, model_weight_oms_flower_5) # ns ttt
  anova(model_weight_oms_flower_3b, model_weight_oms_flower_6) # ns pl
  
  # effect of ratio_q on rs proxies ----
  ## aborted ----
  model_ab_q_flower_0 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ ratio_q * ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_1 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ ratio_q + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_2 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ ratio_q * ttt + pl + ttt + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_3 <- lme4::glmer(data = data_final_flower, cbind(nb_seed,nb_ab) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_3a <- lme4::glmer(data = data_final_flower |> filter (!is.na(ratio_q)), cbind(nb_seed,nb_ab) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_3b <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_seed,nb_ab) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_4 <- lme4::glmer(data = data_final_flower|> filter (!is.na(ratio_q)), cbind(nb_seed,nb_ab) ~ ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_5 <- lme4::glmer(data = data_final_flower , cbind(nb_seed,nb_ab) ~ ratio_q + pl + (1|session) + (1|session:id), family = "binomial")
  model_ab_q_flower_6 <- lme4::glmer(data = data_final_flower  |> filter (!is.na(pl)), cbind(nb_seed,nb_ab) ~ ratio_q + ttt + (1|session) + (1|session:id), family = "binomial")
  anova(model_ab_q_flower_0, model_ab_q_flower_1) # ns ratio_q * ttt
  anova(model_ab_q_flower_0, model_ab_q_flower_2) # ns pl * ttt
  anova(model_ab_q_flower_3a, model_ab_q_flower_4) # ns ratio_q
  anova(model_ab_q_flower_3, model_ab_q_flower_5) # ns sign ttt
  anova(model_ab_q_flower_3b, model_ab_q_flower_6) # ns pl

  ## seed-set ----
  model_ss_q_flower_0 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q * ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_1 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_2 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q * ttt + pl + ttt + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_3 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_3a <- lme4::glmer(data = data_final_flower |> filter (!is.na(ratio_q)), cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_3b <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_4 <- lme4::glmer(data = data_final_flower |> filter (!is.na(ratio_q)), cbind(nb_seed, nb_ov - nb_seed) ~ ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_5 <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q + pl + (1|session) + (1|session:id), family = "binomial")
  model_ss_q_flower_6 <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q + ttt + (1|session) + (1|session:id), family = "binomial")
  anova(model_ss_q_flower_0, model_ss_q_flower_1) # ns ratio_q * ttt
  anova(model_ss_q_flower_0, model_ss_q_flower_2) # marg sign pl * ttt
  anova(model_ss_q_flower_3a, model_ss_q_flower_4) # ns ratio_q
  anova(model_ss_q_flower_3, model_ss_q_flower_5) # ns ttt
  anova(model_ss_q_flower_3b, model_ss_q_flower_6) # sign pl

  data_final_flower <- data_final_flower |>
    mutate(ttt=forcats::fct_relevel(ttt, "high"))
  model_ss_q_flower_best <- lme4::glmer(data = data_final_flower, cbind(nb_seed, nb_ov - nb_seed) ~ ratio_q + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  summary(model_ss_q_flower_best)
  
  ## germ ----
  model_germ_q_flower_0 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q * ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_1 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q + ttt + pl * ttt + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_2 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q * ttt + pl + ttt + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_3 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_3a <- lme4::glmer(data = data_final_flower |> filter (!is.na(ratio_q)), cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_3b <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q + ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_4 <- lme4::glmer(data = data_final_flower |> filter (!is.na(ratio_q)), cbind(nb_germ, nb_sown - nb_germ) ~ ttt + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_5 <- lme4::glmer(data = data_final_flower, cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q + pl + (1|session) + (1|session:id), family = "binomial")
  model_germ_q_flower_6 <- lme4::glmer(data = data_final_flower |> filter (!is.na(pl)), cbind(nb_germ, nb_sown - nb_germ) ~ ratio_q + ttt + (1|session) + (1|session:id), family = "binomial")
  anova(model_germ_q_flower_0, model_germ_q_flower_1) # ns ratio_q * ttt
  anova(model_germ_q_flower_0, model_germ_q_flower_2) # ns pl * ttt
  anova(model_germ_q_flower_3a, model_germ_q_flower_4) # ns ratio_q
  anova(model_germ_q_flower_3, model_germ_q_flower_5) # ns ttt
  anova(model_germ_q_flower_3b, model_germ_q_flower_6) # ns pl

  ## weight ----
  model_weight_q_flower_0 <- lme4::lmer(data = data_final_flower, seed_weight ~ ratio_q * ttt + pl * ttt + (1|session) + (1|session:id))
  model_weight_q_flower_1 <- lme4::lmer(data = data_final_flower, seed_weight ~ ratio_q + ttt + pl * ttt + (1|session) + (1|session:id))
  model_weight_q_flower_2 <- lme4::lmer(data = data_final_flower, seed_weight ~ ratio_q * ttt + pl + ttt + (1|session) + (1|session:id))
  model_weight_q_flower_3 <- lme4::lmer(data = data_final_flower, seed_weight ~ ratio_q + ttt + pl + (1|session) + (1|session:id))
  model_weight_q_flower_3a <- lme4::lmer(data = data_final_flower |> filter (!is.na(ratio_q)), seed_weight ~ ratio_q + ttt + pl + (1|session) + (1|session:id))
  model_weight_q_flower_3b <- lme4::lmer(data = data_final_flower |> filter (!is.na(pl)), seed_weight ~ ratio_q + ttt + pl + (1|session) + (1|session:id))
  model_weight_q_flower_4 <- lme4::lmer(data = data_final_flower |> filter (!is.na(ratio_q)), seed_weight ~ ttt + pl + (1|session) + (1|session:id))
  model_weight_q_flower_5 <- lme4::lmer(data = data_final_flower, seed_weight ~ ratio_q + pl + (1|session) + (1|session:id))
  model_weight_q_flower_6 <- lme4::lmer(data = data_final_flower |> filter (!is.na(pl)), seed_weight ~ ratio_q + ttt + (1|session) + (1|session:id))
  anova(model_weight_q_flower_0, model_weight_q_flower_1) # ns ratio_q * ttt
  anova(model_weight_q_flower_0, model_weight_q_flower_2) # sign pl * ttt
  anova(model_weight_q_flower_3a, model_weight_q_flower_4) # ns ratio_q
  anova(model_weight_q_flower_3, model_weight_q_flower_5) # ns ttt
  anova(model_weight_q_flower_3b, model_weight_q_flower_6) # ns pl

  data_final_flower <- data_final_flower |>
    mutate(ttt=forcats::fct_relevel(ttt, "high"))
  model_weight_q_flower_best <- lmerTest::lmer(data = data_final_flower, seed_weight ~ ratio_q + ttt + pl * ttt + (1|session) + (1|session:id))  
  summary(model_weight_q_flower_best)
  
  return()
}

#' Statistic models: id scale
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

get_stat_id <- function(data_final_id, include_contact = FALSE){

  
  if(include_contact){
    
    method_name <- "w_contact"
    
    # with contact ----
    ## effect of oms on q ratio ----
    model_q_oms_id_0 <- lme4::lmer(data = data_final_id, ratio_q ~ oms * ttt + contact * ttt + (1|session))
    model_q_oms_id_1 <- lme4::lmer(data = data_final_id, ratio_q ~ oms + ttt + contact * ttt + (1|session))
    model_q_oms_id_2 <- lme4::lmer(data = data_final_id, ratio_q ~ oms * ttt + contact + ttt + (1|session))
    model_q_oms_id_3 <- lme4::lmer(data = data_final_id, ratio_q ~ oms + ttt + contact + (1|session))
    model_q_oms_id_4 <- lme4::lmer(data = data_final_id, ratio_q ~ ttt + contact + (1|session))
    model_q_oms_id_5 <- lme4::lmer(data = data_final_id, ratio_q ~ oms + contact + (1|session))
    model_q_oms_id_6 <- lme4::lmer(data = data_final_id, ratio_q ~ oms + ttt + (1|session))
    lrt_interoms_q_oms <- anova(model_q_oms_id_0, model_q_oms_id_1) |> as.data.frame() # oms * ttt
    lrt_intercontact_q_oms <- anova(model_q_oms_id_0, model_q_oms_id_1) |> as.data.frame() # contact * ttt
    lrt_oms_q_oms <- anova(model_q_oms_id_3, model_q_oms_id_4) |> as.data.frame() # oms
    lrt_ttt_q_oms <-anova(model_q_oms_id_3, model_q_oms_id_5) |> as.data.frame() # ttt
    lrt_contact_q_oms <-anova(model_q_oms_id_3, model_q_oms_id_6) |> as.data.frame() # contact
    
    lrt_q_oms <- lrt_interoms_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_intercontact_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_oms_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |>      mutate(method = method_name,
             proxy = "q",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_q_oms_id_best <- lmerTest::lmer(data = data_final_id, ratio_q ~ oms * ttt + contact * ttt + (1|session))
    estimate_low_q_oms <- summary(model_q_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_q_oms_id_best <- lmerTest::lmer(data = data_final_id, ratio_q ~ oms * ttt + contact * ttt + (1|session))
    estimate_medium_q_oms <- summary(model_q_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_q_oms_id_best <- lmerTest::lmer(data = data_final_id, ratio_q ~ oms * ttt + contact * ttt + (1|session))
    estimate_high_q_oms <- summary(model_q_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_q_oms <- estimate_low_q_oms |>
      bind_rows(estimate_medium_q_oms) |>
      bind_rows(estimate_high_q_oms) |>
      rename(pvalue = "Pr(>|t|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ## effect of oms on rs proxies ----
    ### aborted ----
    model_ab_oms_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    model_ab_oms_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms + ttt + contact * ttt + (1|session), family = "binomial")
    model_ab_oms_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + contact + ttt + (1|session), family = "binomial")
    model_ab_oms_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms + ttt + contact + (1|session), family = "binomial")
    model_ab_oms_id_4 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ttt + contact + (1|session), family = "binomial")
    model_ab_oms_id_5 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms + contact + (1|session), family = "binomial")
    model_ab_oms_id_6 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms + ttt + (1|session), family = "binomial")
    lrt_interoms_ab_oms <- anova(model_ab_oms_id_0, model_ab_oms_id_1) |> as.data.frame() # oms * ttt
    lrt_intercontact_ab_oms <- anova(model_ab_oms_id_0, model_ab_oms_id_1) |> as.data.frame() # contact * ttt
    lrt_oms_ab_oms <- anova(model_ab_oms_id_3, model_ab_oms_id_4) |> as.data.frame() # oms
    lrt_ttt_ab_oms <-anova(model_ab_oms_id_3, model_ab_oms_id_5) |> as.data.frame() # ttt
    lrt_contact_ab_oms <-anova(model_ab_oms_id_3, model_ab_oms_id_6) |> as.data.frame() # contact
    
    lrt_ab_oms <- lrt_interoms_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_intercontact_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_oms_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |>      
      mutate(method = method_name,
             proxy = "ab",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ab_oms_id_best <-  lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_low_ab_oms <- summary(model_ab_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ab_oms_id_best <-  lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_medium_ab_oms <- summary(model_ab_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ab_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_high_ab_oms <- summary(model_ab_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_ab_oms <- estimate_low_ab_oms |>
      bind_rows(estimate_medium_ab_oms) |>
      bind_rows(estimate_high_ab_oms) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ### seed-set ----
    model_ss_oms_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    model_ss_oms_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms + ttt + contact * ttt + (1|session), family = "binomial")
    model_ss_oms_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + contact + ttt + (1|session), family = "binomial")
    model_ss_oms_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms + ttt + contact + (1|session), family = "binomial")
    model_ss_oms_id_4 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ttt + contact + (1|session), family = "binomial")
    model_ss_oms_id_5 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms + contact + (1|session), family = "binomial")
    model_ss_oms_id_6 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms + ttt + (1|session), family = "binomial")
    lrt_interoms_ss_oms <- anova(model_ss_oms_id_0, model_ss_oms_id_1) |> as.data.frame() # oms * ttt
    lrt_intercontact_ss_oms <- anova(model_ss_oms_id_0, model_ss_oms_id_1) |> as.data.frame() # contact * ttt
    lrt_oms_ss_oms <- anova(model_ss_oms_id_3, model_ss_oms_id_4) |> as.data.frame() # oms
    lrt_ttt_ss_oms <-anova(model_ss_oms_id_3, model_ss_oms_id_5) |> as.data.frame() # ttt
    lrt_contact_ss_oms <-anova(model_ss_oms_id_3, model_ss_oms_id_6) |> as.data.frame() # contact
    
    lrt_ss_oms <- lrt_interoms_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_intercontact_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_oms_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |>
      mutate(method = method_name,
             proxy = "ss",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ss_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_low_ss_oms <- summary(model_ss_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ss_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_medium_ss_oms <- summary(model_ss_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ss_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_high_ss_oms <- summary(model_ss_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_ss_oms <- estimate_low_ss_oms |>
      bind_rows(estimate_medium_ss_oms) |>
      bind_rows(estimate_high_ss_oms) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ### germ ----
    model_germ_oms_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    model_germ_oms_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms + ttt + contact * ttt + (1|session), family = "binomial")
    model_germ_oms_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + contact + ttt + (1|session), family = "binomial")
    model_germ_oms_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms + ttt + contact + (1|session), family = "binomial")
    model_germ_oms_id_4 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ttt + contact + (1|session), family = "binomial")
    model_germ_oms_id_5 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms + contact + (1|session), family = "binomial")
    model_germ_oms_id_6 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms + ttt + (1|session), family = "binomial")
    lrt_interoms_germ_oms <- anova(model_germ_oms_id_0, model_germ_oms_id_1) |> as.data.frame() # oms * ttt
    lrt_intercontact_germ_oms <- anova(model_germ_oms_id_0, model_germ_oms_id_1) |> as.data.frame() # contact * ttt
    lrt_oms_germ_oms <- anova(model_germ_oms_id_3, model_germ_oms_id_4) |> as.data.frame() # oms
    lrt_ttt_germ_oms <-anova(model_germ_oms_id_3, model_germ_oms_id_5) |> as.data.frame() # ttt
    lrt_contact_germ_oms <-anova(model_germ_oms_id_3, model_germ_oms_id_6) |> as.data.frame() # contact
    
    lrt_germ_oms <- lrt_interoms_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_intercontact_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_oms_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |> 
      mutate(method = method_name,proxy = "germ",
             response = "oms",
             .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_germ_oms_id_best <-lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_low_germ_oms <- summary(model_germ_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_germ_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_medium_germ_oms <- summary(model_germ_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_germ_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_high_germ_oms <- summary(model_germ_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_germ_oms <- estimate_low_germ_oms |>
      bind_rows(estimate_medium_germ_oms) |>
      bind_rows(estimate_high_germ_oms) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ### weight ----
    model_weight_oms_id_0 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + contact * ttt + (1|session))
    model_weight_oms_id_1 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms + ttt + contact * ttt + (1|session))
    model_weight_oms_id_2 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + contact + ttt + (1|session))
    model_weight_oms_id_3 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms + ttt + contact + (1|session))
    model_weight_oms_id_4 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ttt + contact + (1|session))
    model_weight_oms_id_5 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms + contact + (1|session))
    model_weight_oms_id_6 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms + ttt + (1|session))
    lrt_interoms_weight_oms <- anova(model_weight_oms_id_0, model_weight_oms_id_1) |> as.data.frame() # oms * ttt
    lrt_intercontact_weight_oms <- anova(model_weight_oms_id_0, model_weight_oms_id_1) |> as.data.frame() # contact * ttt
    lrt_oms_weight_oms <- anova(model_weight_oms_id_3, model_weight_oms_id_4) |> as.data.frame() # oms
    lrt_ttt_weight_oms <-anova(model_weight_oms_id_3, model_weight_oms_id_5) |> as.data.frame() # ttt
    lrt_contact_weight_oms <-anova(model_weight_oms_id_3, model_weight_oms_id_6) |> as.data.frame() # contact
    
    lrt_weight_oms <- lrt_interoms_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_intercontact_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_oms_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |>
      mutate(method = method_name,proxy = "weight",
             response = "oms", 
             .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_weight_oms_id_best <-lmerTest::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + contact * ttt + (1|session))
    estimate_low_weight_oms <- summary(model_weight_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_weight_oms_id_best <- lmerTest::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + contact * ttt + (1|session))
    estimate_medium_weight_oms <- summary(model_weight_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_weight_oms_id_best <- lmerTest::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + contact * ttt + (1|session))
    estimate_high_weight_oms <- summary(model_weight_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("oms","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_weight_oms <- estimate_low_weight_oms |>
      bind_rows(estimate_medium_weight_oms) |>
      bind_rows(estimate_high_weight_oms) |>
      rename(pvalue = "Pr(>|t|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ## effect of ratio_q on rs proxies ----
    ### aborted ----
    model_ab_q_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    model_ab_q_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + ttt + contact * ttt + (1|session), family = "binomial")
    model_ab_q_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + contact + ttt + (1|session), family = "binomial")
    model_ab_q_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + ttt + contact + (1|session), family = "binomial")
    model_ab_q_id_3a <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + ttt + contact + (1|session), family = "binomial")
    model_ab_q_id_4 <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ttt + contact + (1|session), family = "binomial")
    model_ab_q_id_5 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + contact + (1|session), family = "binomial")
    model_ab_q_id_6 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + ttt + (1|session), family = "binomial")
    lrt_interq_ab_q <- anova(model_ab_q_id_0, model_ab_q_id_1) |> as.data.frame() # ratio_q * ttt
    lrt_intercontact_ab_q <- anova(model_ab_q_id_0, model_ab_q_id_1) |> as.data.frame() # contact * ttt
    lrt_q_ab_q <- anova(model_ab_q_id_3a, model_ab_q_id_4) |> as.data.frame() # ratio_q
    lrt_ttt_ab_q <-anova(model_ab_q_id_3, model_ab_q_id_5) |> as.data.frame() # ttt
    lrt_contact_ab_q <-anova(model_ab_q_id_3, model_ab_q_id_6) |> as.data.frame() # contact
    
    lrt_ab_q <- lrt_interq_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_intercontact_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_q_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |>      
      mutate(method = method_name,
             proxy = "ab",
             response = "ratio_q", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ab_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_low_ab_q <- summary(model_ab_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ab_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_medium_ab_q <- summary(model_ab_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ab_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_high_ab_q <- summary(model_ab_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_ab_q <- estimate_low_ab_q |>
      bind_rows(estimate_medium_ab_q) |>
      bind_rows(estimate_high_ab_q) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ### seed-set ----
    model_ss_q_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    model_ss_q_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + ttt + contact * ttt + (1|session), family = "binomial")
    model_ss_q_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + contact + ttt + (1|session), family = "binomial")
    model_ss_q_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + ttt + contact + (1|session), family = "binomial")
    model_ss_q_id_3a <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + ttt + contact + (1|session), family = "binomial")
    model_ss_q_id_4 <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ttt + contact + (1|session), family = "binomial")
    model_ss_q_id_5 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + contact + (1|session), family = "binomial")
    model_ss_q_id_6 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + ttt + (1|session), family = "binomial")
    lrt_interq_ss_q <- anova(model_ss_q_id_0, model_ss_q_id_1) |> as.data.frame() # ratio_q * ttt
    lrt_intercontact_ss_q <- anova(model_ss_q_id_0, model_ss_q_id_1) |> as.data.frame() # contact * ttt
    lrt_q_ss_q <- anova(model_ss_q_id_3a, model_ss_q_id_4) |> as.data.frame() # ratio_q
    lrt_ttt_ss_q <-anova(model_ss_q_id_3, model_ss_q_id_5) |> as.data.frame() # ttt
    lrt_contact_ss_q <-anova(model_ss_q_id_3, model_ss_q_id_6) |> as.data.frame() # contact
    
    lrt_ss_q <- lrt_interq_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_intercontact_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_q_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |>
      mutate(method = method_name,
             proxy = "ss",
             response = "ratio_q", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ss_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_low_ss_q <- summary(model_ss_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ss_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_medium_ss_q <- summary(model_ss_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ss_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_high_ss_q <- summary(model_ss_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_ss_q <- estimate_low_ss_q |>
      bind_rows(estimate_medium_ss_q) |>
      bind_rows(estimate_high_ss_q) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ### germ ----
    model_germ_q_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    model_germ_q_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + ttt + contact * ttt + (1|session), family = "binomial")
    model_germ_q_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + contact + ttt + (1|session), family = "binomial")
    model_germ_q_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + ttt + contact + (1|session), family = "binomial")
    model_germ_q_id_3a <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + ttt + contact + (1|session), family = "binomial")
    model_germ_q_id_4 <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ttt + contact + (1|session), family = "binomial")
    model_germ_q_id_5 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + contact + (1|session), family = "binomial")
    model_germ_q_id_6 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + ttt + (1|session), family = "binomial")
    lrt_interq_germ_q <- anova(model_germ_q_id_0, model_germ_q_id_1) |> as.data.frame() # ratio_q * ttt
    lrt_intercontact_germ_q <- anova(model_germ_q_id_0, model_germ_q_id_1) |> as.data.frame() # contact * ttt
    lrt_q_germ_q <- anova(model_germ_q_id_3a, model_germ_q_id_4) |> as.data.frame() # ratio_q
    lrt_ttt_germ_q <-anova(model_germ_q_id_3, model_germ_q_id_5) |> as.data.frame() # ttt
    lrt_contact_germ_q <-anova(model_germ_q_id_3, model_germ_q_id_6) |> as.data.frame() # contact
    
    lrt_germ_q <- lrt_interq_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_intercontact_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_q_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |> 
      mutate(method = method_name,proxy = "germ",
             response = "ratio_q",
             .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_germ_q_id_best <-lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_low_germ_q <- summary(model_germ_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_germ_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_medium_germ_q <- summary(model_germ_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_germ_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + contact * ttt + (1|session), family = "binomial")
    estimate_high_germ_q <- summary(model_germ_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_germ_q <- estimate_low_germ_q |>
      bind_rows(estimate_medium_germ_q) |>
      bind_rows(estimate_high_germ_q) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    ### weight ----
    model_weight_q_id_0 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + contact * ttt + (1|session))
    model_weight_q_id_1 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q + ttt + contact * ttt + (1|session))
    model_weight_q_id_2 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + contact + ttt + (1|session))
    model_weight_q_id_3 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q + ttt + contact + (1|session))
    model_weight_q_id_3a <- lme4::lmer(data = data_final_id |> filter(!(is.na(ratio_q))), mean_seed_weight ~ ratio_q + ttt + contact + (1|session))
    model_weight_q_id_4 <- lme4::lmer(data = data_final_id |> filter(!(is.na(ratio_q))), mean_seed_weight ~ ttt + contact + (1|session))
    model_weight_q_id_5 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q + contact + (1|session))
    model_weight_q_id_6 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q + ttt + (1|session))
    lrt_interq_weight_q <- anova(model_weight_q_id_0, model_weight_q_id_1) |> as.data.frame() # ratio_q * ttt
    lrt_intercontact_weight_q <- anova(model_weight_q_id_0, model_weight_q_id_1) |> as.data.frame() # contact * ttt
    lrt_q_weight_q <- anova(model_weight_q_id_3a, model_weight_q_id_4) |> as.data.frame() # ratio_q
    lrt_ttt_weight_q <-anova(model_weight_q_id_3, model_weight_q_id_5) |> as.data.frame() # ttt
    lrt_contact_weight_q <-anova(model_weight_q_id_3, model_weight_q_id_6) |> as.data.frame() # contact
    
    lrt_weight_q <- lrt_interq_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_intercontact_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "intercontact", .before = 1)) |>
      bind_rows(lrt_q_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      bind_rows(lrt_contact_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "contact", .before = 1)) |>
      mutate(method = method_name,proxy = "weight",
             response = "ratio_q", 
             .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_weight_q_id_best <-lmerTest::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + contact * ttt + (1|session))
    estimate_low_weight_q <- summary(model_weight_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_weight_q_id_best <- lmerTest::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + contact * ttt + (1|session))
    estimate_medium_weight_q <- summary(model_weight_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_weight_q_id_best <- lmerTest::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + contact * ttt + (1|session))
    estimate_high_weight_q <- summary(model_weight_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname %in% c("ratio_q","contact")) |> mutate(ttt = "high", .before = 1)
    
    estimate_weight_q <- estimate_low_weight_q |>
      bind_rows(estimate_medium_weight_q) |>
      bind_rows(estimate_high_weight_q) |>
      rename(pvalue = "Pr(>|t|)") |> 
      mutate(method = method_name,
             proxy = "q", .before = 1)
    
    lrt_table <- lrt_q_oms |>
      bind_rows(lrt_ab_oms) |>
      bind_rows(lrt_ss_oms) |>
      bind_rows(lrt_germ_oms) |>
      bind_rows(lrt_weight_oms) |>
      bind_rows(lrt_ab_q) |>
      bind_rows(lrt_ss_q) |>
      bind_rows(lrt_germ_q) |>
      bind_rows(lrt_weight_q) |>
      rename(chisq = Chisq,
             df = Df,
             pvalue = "Pr(>Chisq)") |>
      mutate(level = "id", .before = 1)
    
    estimate_table <- estimate_q_oms |>
      bind_rows(estimate_ab_oms) |>
      bind_rows(estimate_ss_oms) |>
      bind_rows(estimate_germ_oms) |>
      bind_rows(estimate_weight_oms) |>
      bind_rows(estimate_ab_q) |>
      bind_rows(estimate_ss_q) |>
      bind_rows(estimate_germ_q) |>
      bind_rows(estimate_weight_q) |>
      select(-c("t value","z value")) |>
      rename(estimate = Estimate,
             se = "Std. Error") |>
      mutate(level = "id", .before = 1)
    
    
  }else{
    
    method_name <- "wo_contact"
    
    # without contact ----
    ## effect of oms on q ratio ----
    model_q_oms_id_0 <- lme4::lmer(data = data_final_id, ratio_q ~ oms * ttt + (1|session))
    model_q_oms_id_1 <- lme4::lmer(data = data_final_id, ratio_q ~ oms + ttt + (1|session))
    model_q_oms_id_2 <- lme4::lmer(data = data_final_id, ratio_q ~ ttt + (1|session))
    model_q_oms_id_3 <- lme4::lmer(data = data_final_id, ratio_q ~ oms + (1|session))
    lrt_interoms_q_oms <- anova(model_q_oms_id_0, model_q_oms_id_1) |> as.data.frame() # ns oms * ttt
    lrt_oms_q_oms <- anova(model_q_oms_id_1, model_q_oms_id_2) |> as.data.frame() # marg sin oms
    lrt_ttt_q_oms <-anova(model_q_oms_id_1, model_q_oms_id_3) |> as.data.frame() # ns ttt
      
    lrt_q_oms <- lrt_interoms_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_oms_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_q_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "q",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_q_oms_id_best <- lmerTest::lmer(data = data_final_id, ratio_q ~ oms * ttt + (1|session))
    estimate_low_q_oms <- summary(model_q_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_q_oms_id_best <- lmerTest::lmer(data = data_final_id, ratio_q ~ oms * ttt + (1|session))
    estimate_medium_q_oms <- summary(model_q_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_q_oms_id_best <- lmerTest::lmer(data = data_final_id, ratio_q ~ oms * ttt + (1|session))
    estimate_high_q_oms <- summary(model_q_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "high", .before = 1)
    
    estimate_q_oms <- estimate_low_q_oms |>
      bind_rows(estimate_medium_q_oms) |>
      bind_rows(estimate_high_q_oms) |>
      rename(pvalue = "Pr(>|t|)") |> 
      mutate(method = method_name,
             proxy = "q",
             response = "oms", .before = 1)
    
    ## effect of oms on rs proxies ----
    ### aborted ----
    model_ab_oms_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + (1|session), family = "binomial")
    model_ab_oms_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms + ttt + (1|session), family = "binomial")
    model_ab_oms_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ttt + (1|session), family = "binomial")
    model_ab_oms_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms + (1|session), family = "binomial")
    lrt_interoms_ab_oms <- anova(model_ab_oms_id_0, model_ab_oms_id_1) |> as.data.frame() # ns oms * ttt
    lrt_oms_ab_oms <- anova(model_ab_oms_id_1, model_ab_oms_id_2) |> as.data.frame() # marg sin oms
    lrt_ttt_ab_oms <-anova(model_ab_oms_id_1, model_ab_oms_id_3) |> as.data.frame() # ns ttt
    
    lrt_ab_oms <- lrt_interoms_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_oms_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_ab_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "ab",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ab_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + (1|session), family = "binomial")
    estimate_low_ab_oms <- summary(model_ab_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ab_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + (1|session), family = "binomial")
    estimate_medium_ab_oms <- summary(model_ab_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ab_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ oms * ttt + (1|session), family = "binomial")
    estimate_high_ab_oms <- summary(model_ab_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "high", .before = 1)
    
    estimate_ab_oms <- estimate_low_ab_oms |>
      bind_rows(estimate_medium_ab_oms) |>
      bind_rows(estimate_high_ab_oms) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "ab",
             response = "oms", .before = 1)
    
    ### seed-set ----
    model_ss_oms_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + (1|session), family = "binomial")
    model_ss_oms_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms + ttt + (1|session), family = "binomial")
    model_ss_oms_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ttt + (1|session), family = "binomial")
    model_ss_oms_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms + (1|session), family = "binomial")
    lrt_interoms_ss_oms <- anova(model_ss_oms_id_0, model_ss_oms_id_1) |> as.data.frame() # ns oms * ttt
    lrt_oms_ss_oms <- anova(model_ss_oms_id_1, model_ss_oms_id_2) |> as.data.frame() # marg sin oms
    lrt_ttt_ss_oms <-anova(model_ss_oms_id_1, model_ss_oms_id_3) |> as.data.frame() # ns ttt
    
    lrt_ss_oms <- lrt_interoms_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_oms_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_ss_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "ss",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ss_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + (1|session), family = "binomial")
    estimate_low_ss_oms <- summary(model_ss_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ss_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + (1|session), family = "binomial")
    estimate_medium_ss_oms <- summary(model_ss_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ss_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ oms * ttt + (1|session), family = "binomial")
    estimate_high_ss_oms <- summary(model_ss_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "high", .before = 1)
    
    estimate_ss_oms <- estimate_low_ss_oms |>
      bind_rows(estimate_medium_ss_oms) |>
      bind_rows(estimate_high_ss_oms) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "ss",
             response = "oms", .before = 1)
    
    ### germ ----
    model_germ_oms_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + (1|session), family = "binomial")
    model_germ_oms_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms + ttt + (1|session), family = "binomial")
    model_germ_oms_id_2 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ttt + (1|session), family = "binomial")
    model_germ_oms_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms + (1|session), family = "binomial")
    lrt_interoms_germ_oms <- anova(model_germ_oms_id_0, model_germ_oms_id_1) |> as.data.frame() # ns oms * ttt
    lrt_oms_germ_oms <- anova(model_germ_oms_id_1, model_germ_oms_id_2) |> as.data.frame() # marg sin oms
    lrt_ttt_germ_oms <-anova(model_germ_oms_id_1, model_germ_oms_id_3) |> as.data.frame() # ns ttt
    
    lrt_germ_oms <- lrt_interoms_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_oms_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_germ_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "germ",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_germ_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + (1|session), family = "binomial")
    estimate_low_germ_oms <- summary(model_germ_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_germ_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + (1|session), family = "binomial")
    estimate_medium_germ_oms <- summary(model_germ_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_germ_oms_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ oms * ttt + (1|session), family = "binomial")
    estimate_high_germ_oms <- summary(model_germ_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "high", .before = 1)
    
    estimate_germ_oms <- estimate_low_germ_oms |>
      bind_rows(estimate_medium_germ_oms) |>
      bind_rows(estimate_high_germ_oms) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "germ",
             response = "oms", .before = 1)
    
    ### weight ----
    model_weight_oms_id_0 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + (1|session))
    model_weight_oms_id_1 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms + ttt + (1|session))
    model_weight_oms_id_2 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ttt + (1|session))
    model_weight_oms_id_3 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ oms + (1|session))
    lrt_interoms_weight_oms <- anova(model_weight_oms_id_0, model_weight_oms_id_1) |> as.data.frame() # ns oms * ttt
    lrt_oms_weight_oms <- anova(model_weight_oms_id_1, model_weight_oms_id_2) |> as.data.frame() # marg sin oms
    lrt_ttt_weight_oms <-anova(model_weight_oms_id_1, model_weight_oms_id_3) |> as.data.frame() # ns ttt
    
    lrt_weight_oms <- lrt_interoms_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interoms", .before = 1) |>
      bind_rows(lrt_oms_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "oms", .before = 1)) |>
      bind_rows(lrt_ttt_weight_oms |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "weight",
             response = "oms", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_weight_oms_id_best <-  lmerTest::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + (1|session))
    estimate_low_weight_oms <- summary(model_weight_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_weight_oms_id_best <-  lmerTest::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + (1|session))
    estimate_medium_weight_oms <- summary(model_weight_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_weight_oms_id_best <-  lmerTest::lmer(data = data_final_id, mean_seed_weight ~ oms * ttt + (1|session))
    estimate_high_weight_oms <- summary(model_weight_oms_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "oms") |> mutate(ttt = "high", .before = 1)
    
    estimate_weight_oms <- estimate_low_weight_oms |>
      bind_rows(estimate_medium_weight_oms) |>
      bind_rows(estimate_high_weight_oms) |>
      rename(pvalue = "Pr(>|t|)") |> 
      mutate(method = method_name,
             proxy = "weight",
             response = "oms", .before = 1)
    
    ## effect of ratio_q on rs proxies ----
    ### aborted ----
    model_ab_q_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + (1|session), family = "binomial")
    model_ab_q_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + ttt + (1|session), family = "binomial")
    model_ab_q_id_1b <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)) , cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + ttt + (1|session), family = "binomial")
    model_ab_q_id_2 <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)) , cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ttt + (1|session), family = "binomial")
    model_ab_q_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q + (1|session), family = "binomial")
    lrt_interq_ab_q <- anova(model_ab_q_id_0, model_ab_q_id_1) |> as.data.frame() # ns ratio_q * ttt
    lrt_q_ab_q <- anova(model_ab_q_id_1b, model_ab_q_id_2) |> as.data.frame() # marg sin ratio_q
    lrt_ttt_ab_q <-anova(model_ab_q_id_1, model_ab_q_id_3) |> as.data.frame() # ns ttt
    
    lrt_ab_q <- lrt_interq_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_q_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_ab_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "ab",
             response = "ratio_q", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ab_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_low_ab_q <- summary(model_ab_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ab_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_medium_ab_q <- summary(model_ab_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ab_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ab,nb_ab_sum_ab) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_high_ab_q <- summary(model_ab_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "high", .before = 1)
    
    estimate_ab_q <- estimate_low_ab_q |>
      bind_rows(estimate_medium_ab_q) |>
      bind_rows(estimate_high_ab_q) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "ab",
             response = "ratio_q", .before = 1)
    
    ### seed-set ----
    model_ss_q_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + (1|session), family = "binomial")
    model_ss_q_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + ttt + (1|session), family = "binomial")
    model_ss_q_id_1b <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + ttt + (1|session), family = "binomial")
    model_ss_q_id_2 <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ttt + (1|session), family = "binomial")
    model_ss_q_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q + (1|session), family = "binomial")
    lrt_interq_ss_q <- anova(model_ss_q_id_0, model_ss_q_id_1) |> as.data.frame() # ns ratio_q * ttt
    lrt_q_ss_q <- anova(model_ss_q_id_1b, model_ss_q_id_2) |> as.data.frame() # marg sin ratio_q
    lrt_ttt_ss_q <-anova(model_ss_q_id_1, model_ss_q_id_3) |> as.data.frame() # ns ttt
    
    lrt_ss_q <- lrt_interq_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_q_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_ss_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "ss",
             response = "ratio_q", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_ss_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_low_ss_q <- summary(model_ss_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_ss_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_medium_ss_q <- summary(model_ss_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_ss_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_seed_sum_ss, nb_ov_sum_ss - nb_seed_sum_ss) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_high_ss_q <- summary(model_ss_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "high", .before = 1)
    
    estimate_ss_q <- estimate_low_ss_q |>
      bind_rows(estimate_medium_ss_q) |>
      bind_rows(estimate_high_ss_q) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "ss",
             response = "ratio_q", .before = 1)
    
    ### germ ----
    model_germ_q_id_0 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + (1|session), family = "binomial")
    model_germ_q_id_1 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + ttt + (1|session), family = "binomial")
    model_germ_q_id_1b <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + ttt + (1|session), family = "binomial")
    model_germ_q_id_2 <- lme4::glmer(data = data_final_id |> filter(!is.na(ratio_q)), cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ttt + (1|session), family = "binomial")
    model_germ_q_id_3 <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q + (1|session), family = "binomial")
    lrt_interq_germ_q <- anova(model_germ_q_id_0, model_germ_q_id_1) |> as.data.frame() # ns ratio_q * ttt
    lrt_q_germ_q <- anova(model_germ_q_id_1b, model_germ_q_id_2) |> as.data.frame() # marg sin ratio_q
    lrt_ttt_germ_q <-anova(model_germ_q_id_1, model_germ_q_id_3) |> as.data.frame() # ns ttt
    
    lrt_germ_q <- lrt_interq_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_q_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_germ_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "germ",
             response = "ratio_q", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_germ_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_low_germ_q <- summary(model_germ_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_germ_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_medium_germ_q <- summary(model_germ_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_germ_q_id_best <- lme4::glmer(data = data_final_id, cbind(nb_germ_sum_germ, nb_sown_sum_germ - nb_germ_sum_germ) ~ ratio_q * ttt + (1|session), family = "binomial")
    estimate_high_germ_q <- summary(model_germ_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "high", .before = 1)
    
    estimate_germ_q <- estimate_low_germ_q |>
      bind_rows(estimate_medium_germ_q) |>
      bind_rows(estimate_high_germ_q) |>
      rename(pvalue = "Pr(>|z|)") |> 
      mutate(method = method_name,
             proxy = "germ",
             response = "ratio_q", .before = 1)
    
    ### weight ----
    model_weight_q_id_0 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + (1|session))
    model_weight_q_id_1 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q + ttt + (1|session))
    model_weight_q_id_1b <- lme4::lmer(data = data_final_id |> filter(!is.na(ratio_q)), mean_seed_weight ~ ratio_q + ttt + (1|session))
    model_weight_q_id_2 <- lme4::lmer(data = data_final_id |> filter(!is.na(ratio_q)), mean_seed_weight ~ ttt + (1|session))
    model_weight_q_id_3 <- lme4::lmer(data = data_final_id, mean_seed_weight ~ ratio_q + (1|session))
    lrt_interq_weight_q <- anova(model_weight_q_id_0, model_weight_q_id_1) |> as.data.frame() # ns ratio_q * ttt
    lrt_q_weight_q <- anova(model_weight_q_id_1b, model_weight_q_id_2) |> as.data.frame() # marg sin ratio_q
    lrt_ttt_weight_q <-anova(model_weight_q_id_1, model_weight_q_id_3) |> as.data.frame() # ns ttt
    
    lrt_weight_q <- lrt_interq_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "interq", .before = 1) |>
      bind_rows(lrt_q_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ratio_q", .before = 1)) |>
      bind_rows(lrt_ttt_weight_q |> slice(2) |> select(c("Chisq","Df","Pr(>Chisq)")) |> mutate(effect = "ttt", .before = 1)) |>
      mutate(method = method_name,
             proxy = "weight",
             response = "ratio_q", .before = 1)
    
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "low"))
    model_weight_q_id_best <-  lmerTest::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + (1|session))
    estimate_low_weight_q <- summary(model_weight_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "low", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "medium"))
    model_weight_q_id_best <-  lmerTest::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + (1|session))
    estimate_medium_weight_q <- summary(model_weight_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "medium", .before = 1)
    data_final_id <- data_final_id |>
      mutate(ttt=forcats::fct_relevel(ttt, "high"))
    model_weight_q_id_best <-  lmerTest::lmer(data = data_final_id, mean_seed_weight ~ ratio_q * ttt + (1|session))
    estimate_high_weight_q <- summary(model_weight_q_id_best)$coefficients |> as.data.frame() |> rownames_to_column() |>
      filter(rowname == "ratio_q") |> mutate(ttt = "high", .before = 1)
    
    estimate_weight_q <- estimate_low_weight_q |>
      bind_rows(estimate_medium_weight_q) |>
      bind_rows(estimate_high_weight_q) |>
      rename(pvalue = "Pr(>|t|)") |> 
      mutate(method = method_name,
             proxy = "weight",
             response = "ratio_q", .before = 1)
    
    lrt_table <- lrt_q_oms |>
      bind_rows(lrt_ab_oms) |>
      bind_rows(lrt_ss_oms) |>
      bind_rows(lrt_germ_oms) |>
      bind_rows(lrt_weight_oms) |>
      bind_rows(lrt_ab_q) |>
      bind_rows(lrt_ss_q) |>
      bind_rows(lrt_germ_q) |>
      bind_rows(lrt_weight_q) |>
      rename(chisq = Chisq,
             df = Df,
             pvalue = "Pr(>Chisq)")
    
    estimate_table <- estimate_q_oms |>
      bind_rows(estimate_ab_oms) |>
      bind_rows(estimate_ss_oms) |>
      bind_rows(estimate_germ_oms) |>
      bind_rows(estimate_weight_oms) |>
      bind_rows(estimate_ab_q) |>
      bind_rows(estimate_ss_q) |>
      bind_rows(estimate_germ_q) |>
      bind_rows(estimate_weight_q) |>
      select(-c("t value","z value")) |>
      rename(estimate = Estimate,
             se = "Std. Error")
    
  }
  return(list(lrt_table = lrt_table,
              estimate_table = estimate_table))
}

#' #' Read data at the individual level
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export
#' 
#' load_data_id <- function(file_path){
#'   
#'   data_id <- read.table(file_path,head=T) %>%
#'     mutate(poll_treat_factor=as.factor(case_when(poll_treat==1~"low",
#'                                                  poll_treat==2~"medium",
#'                                                  TRUE~"high"))) %>%
#'     mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high"))) %>%
#'     mutate(prop_self=SR_self/(SR_self+SR_fem_out+SR_mal_out),
#'            gam_ov_proxy=nb_flo*nbOv_mean,
#'            SR_out=SR_fem_out+SR_mal_out) %>%
#'     mutate(mean_nb_visit_per_flower=nb_visit/nb_dist_vis)
#'   
#'   return(data_id)
#' }
#' 
#' 
#' #' Get common theme for ggplot
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export
#' get_common_theme <- function(){
#'   
#'   common_theme <- theme(axis.text.x = element_text(size=10,angle=45,hjust=1),
#'                         axis.text.y = element_text(size=8,margin = margin(t = 0, r = 5, b = 0, l = 5)),
#'                         legend.position = "none",
#'                         plot.background = element_blank(),
#'                         panel.grid.major=element_blank(),
#'                         plot.margin = margin(t=5, r=5, b=5, l=5),
#'                         strip.background = element_rect(colour="black", 
#'                                                         linewidth=1.5, linetype="solid"),
#'                         panel.background = element_rect(colour="black", fill="white",linewidth=1),
#'                         axis.title.y = element_text(size=12))
#'   return(common_theme)
#' }
#' 
#' #' Get proportion of non-aborted
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export
#' 
#' 
#' get_prop_non_abort <- function(data_proxy_av, common_theme){
#'   
#'   scaleFUN <- function(x){sprintf("%.1f", x)}
#'   
#'   # gms
#'   
#'   data_proxy_av <- data_proxy_av %>%
#'     mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_gMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS*poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_av_1 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS+poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_av_2 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_av_3 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS+(1|session),family=binomial)
#'   anova(mod_gMS_av_0,mod_gMS_av_1) # inter
#'   anova(mod_gMS_av_1,mod_gMS_av_2) # MS
#'   anova(mod_gMS_av_1,mod_gMS_av_3) # ttt
#'   gMS_av_low_est <- summary(mod_gMS_av_0)$coefficients["gMS","Estimate"] # low
#'   gMS_av_low_se <- summary(mod_gMS_av_0)$coefficients["gMS","Std. Error"] # low
#'   gMS_av_low_pval <- summary(mod_gMS_av_0)$coefficients["gMS","Pr(>|z|)"] # low
#'   data_proxy_av <- data_proxy_av %>%
#'     mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_gMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS*poll_treat_factor+(1|session),family=binomial)
#'   gMS_av_med_est <- summary(mod_gMS_av_0)$coefficients["gMS","Estimate"] # med
#'   gMS_av_med_se <- summary(mod_gMS_av_0)$coefficients["gMS","Std. Error"]# med
#'   gMS_av_med_pval <- summary(mod_gMS_av_0)$coefficients["gMS","Pr(>|z|)"] # med
#'   data_proxy_av <- data_proxy_av %>%
#'     mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_gMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS*poll_treat_factor+(1|session),family=binomial)
#'   gMS_av_high_est <- summary(mod_gMS_av_0)$coefficients["gMS","Estimate"] # high
#'   gMS_av_high_se <- summary(mod_gMS_av_0)$coefficients["gMS","Std. Error"] # high
#'   gMS_av_high_pval <- summary(mod_gMS_av_0)$coefficients["gMS","Pr(>|z|)"]# high
#'   
#'   tab_gMS_av <- tibble(ttt=c("low","medium","high"),
#'                        est=c(gMS_av_low_est,gMS_av_med_est,gMS_av_high_est),
#'                        se=c(gMS_av_low_se,gMS_av_med_se,gMS_av_high_se),
#'                        pval=c(gMS_av_low_pval,gMS_av_med_pval,gMS_av_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=forcats::fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_gMS_av_fem <- ggplot(tab_gMS_av, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     xlab("") +
#'     ylab("\nProportion of non-aborted seeds \nin relation to genetic number of mates")  + 
#'     ylim(-0.3,1.1) +
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_gMS_av_fem
#'   ggsave("figures/plot_gMS_av_fem.jpeg",plot_gMS_av_fem,width=4,height=4, device = png,bg="white")
#'   
#'   # oms
#'   
#'   # PROPORTION NON-ABORTED
#'   data_proxy_av <- data_proxy_av %>%
#'     mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_oMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS*poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_av_1 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS+poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_av_2 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_av_3 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS+(1|session),family=binomial)
#'   
#'   anova(mod_oMS_av_0,mod_oMS_av_1) # inter
#'   anova(mod_oMS_av_1,mod_oMS_av_2) # MS
#'   anova(mod_oMS_av_1,mod_oMS_av_3) # ttt
#'   oMS_av_low_est <- summary(mod_oMS_av_0)$coefficients["oMS","Estimate"] # low
#'   oMS_av_low_se <- summary(mod_oMS_av_0)$coefficients["oMS","Std. Error"] # low
#'   oMS_av_low_pval <- summary(mod_oMS_av_0)$coefficients["oMS","Pr(>|z|)"] # low
#'   data_proxy_av <- data_proxy_av %>%
#'     mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_oMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS*poll_treat_factor+(1|session),family=binomial)
#'   oMS_av_med_est <- summary(mod_oMS_av_0)$coefficients["oMS","Estimate"] # med
#'   oMS_av_med_se <- summary(mod_oMS_av_0)$coefficients["oMS","Std. Error"]# med
#'   oMS_av_med_pval <- summary(mod_oMS_av_0)$coefficients["oMS","Pr(>|z|)"] # med
#'   data_proxy_av <- data_proxy_av %>%
#'     mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_oMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS*poll_treat_factor+(1|session),family=binomial)
#'   oMS_av_high_est <- summary(mod_oMS_av_0)$coefficients["oMS","Estimate"] # high
#'   oMS_av_high_se <- summary(mod_oMS_av_0)$coefficients["oMS","Std. Error"] # high
#'   oMS_av_high_pval <- summary(mod_oMS_av_0)$coefficients["oMS","Pr(>|z|)"]# high
#'   
#'   tab_oMS_av <- tibble(ttt=c("low","medium","high"),
#'                        est=c(oMS_av_low_est,oMS_av_med_est,oMS_av_high_est),
#'                        se=c(oMS_av_low_se,oMS_av_med_se,oMS_av_high_se),
#'                        pval=c(oMS_av_low_pval,oMS_av_med_pval,oMS_av_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=forcats::fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_oMS_av_fem <- ggplot(tab_oMS_av, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     xlab("") +
#'     scale_y_continuous(limits=c(-0.6,0.8),labels=scaleFUN) +
#'     ylab("\nProportion of non-aborted seeds \nin relation to observational number of mates")  + 
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_oMS_av_fem
#'   ggsave("figures/plot_oMS_av_fem.jpeg",plot_oMS_av_fem,width=4,height=4, device = png,bg="white")
#'   
#'   
#'   
#'   return(list(plot_gms_av = plot_gMS_av_fem,
#'               plot_oms_av = plot_oMS_av_fem))
#' }
#' 
#' #' Get seed-set
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export
#' 
#' 
#' get_seed_set <- function(data_proxy_ov, common_theme){
#'   
#'   scaleFUN <- function(x){sprintf("%.1f", x)}
#'   
#'   # gms
#'   
#'   data_proxy_ov <- data_proxy_ov %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_gMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_ss_1 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS+poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_ss_2 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_ss_3 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS+(1|session),family=binomial)
#'   anova(mod_gMS_ss_0,mod_gMS_ss_1) # inter
#'   anova(mod_gMS_ss_1,mod_gMS_ss_2) # MS
#'   anova(mod_gMS_ss_1,mod_gMS_ss_3) # ttt
#'   gMS_ss_low_est <- summary(mod_gMS_ss_0)$coefficients["gMS","Estimate"] # low
#'   gMS_ss_low_se <- summary(mod_gMS_ss_0)$coefficients["gMS","Std. Error"] # low
#'   gMS_ss_low_pval <- summary(mod_gMS_ss_0)$coefficients["gMS","Pr(>|z|)"] # low
#'   data_proxy_ov <- data_proxy_ov %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_gMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
#'   gMS_ss_med_est <- summary(mod_gMS_ss_0)$coefficients["gMS","Estimate"] # med
#'   gMS_ss_med_se <- summary(mod_gMS_ss_0)$coefficients["gMS","Std. Error"]# med
#'   gMS_ss_med_pval <- summary(mod_gMS_ss_0)$coefficients["gMS","Pr(>|z|)"] # med
#'   data_proxy_ov <- data_proxy_ov %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_gMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
#'   gMS_ss_high_est <- summary(mod_gMS_ss_0)$coefficients["gMS","Estimate"] # high
#'   gMS_ss_high_se <- summary(mod_gMS_ss_0)$coefficients["gMS","Std. Error"] # high
#'   gMS_ss_high_pval <- summary(mod_gMS_ss_0)$coefficients["gMS","Pr(>|z|)"]# high
#'   
#'   tab_gMS_ss <- tibble(ttt=c("low","medium","high"),
#'                        est=c(gMS_ss_low_est,gMS_ss_med_est,gMS_ss_high_est),
#'                        se=c(gMS_ss_low_se,gMS_ss_med_se,gMS_ss_high_se),
#'                        pval=c(gMS_ss_low_pval,gMS_ss_med_pval,gMS_ss_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_gMS_ss_fem <- ggplot(tab_gMS_ss, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     ylim(0,0.5) +
#'     xlab("") +
#'     ylab("\nSeed-set in relation to\ngenetic number of mates")  + 
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_gMS_ss_fem
#'   ggsave("figures/plot_gMS_ss_fem.jpeg",plot_gMS_ss_fem,width=4,height=4, device = png,bg="white")
#'   
#'   # oms
#'   
#'   data_proxy_ov <- data_proxy_ov %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_oMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_ss_1 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS+poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_ss_2 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_ss_3 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS+(1|session),family=binomial)
#'   
#'   anova(mod_oMS_ss_0,mod_oMS_ss_1) # inter
#'   anova(mod_oMS_ss_1,mod_oMS_ss_2) # MS
#'   anova(mod_oMS_ss_1,mod_oMS_ss_3) # ttt
#'   oMS_ss_low_est <- summary(mod_oMS_ss_0)$coefficients["oMS","Estimate"] # low
#'   oMS_ss_low_se <- summary(mod_oMS_ss_0)$coefficients["oMS","Std. Error"] # low
#'   oMS_ss_low_pval <- summary(mod_oMS_ss_0)$coefficients["oMS","Pr(>|z|)"] # low
#'   data_proxy_ov <- data_proxy_ov %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_oMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
#'   oMS_ss_med_est <- summary(mod_oMS_ss_0)$coefficients["oMS","Estimate"] # med
#'   oMS_ss_med_se <- summary(mod_oMS_ss_0)$coefficients["oMS","Std. Error"]# med
#'   oMS_ss_med_pval <- summary(mod_oMS_ss_0)$coefficients["oMS","Pr(>|z|)"] # med
#'   data_proxy_ov <- data_proxy_ov %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_oMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
#'   oMS_ss_high_est <- summary(mod_oMS_ss_0)$coefficients["oMS","Estimate"] # high
#'   oMS_ss_high_se <- summary(mod_oMS_ss_0)$coefficients["oMS","Std. Error"] # high
#'   oMS_ss_high_pval <- summary(mod_oMS_ss_0)$coefficients["oMS","Pr(>|z|)"]# high
#'   
#'   tab_oMS_ss <- tibble(ttt=c("low","medium","high"),
#'                        est=c(oMS_ss_low_est,oMS_ss_med_est,oMS_ss_high_est),
#'                        se=c(oMS_ss_low_se,oMS_ss_med_se,oMS_ss_high_se),
#'                        pval=c(oMS_ss_low_pval,oMS_ss_med_pval,oMS_ss_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_oMS_ss_fem <- ggplot(tab_oMS_ss, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     xlab("") +
#'     scale_y_continuous(limits=c(-0.15,0.4),labels=scaleFUN) +
#'     ylab("\nSeed-set in relation to\nobservational number of mates")  + 
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_oMS_ss_fem
#'   ggsave("figures/plot_oMS_ss_fem.jpeg",plot_oMS_ss_fem,width=4,height=4, device = png,bg="white")
#'   
#'   
#'   return(list(plot_gms_ss = plot_gMS_ss_fem,
#'               plot_oms_ss = plot_oMS_ss_fem))
#' }
#' 
#' #' Get proportion germinated seeds
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export
#' 
#' 
#' get_prop_germ <- function(data_id_fem, common_theme){
#'   
#'   scaleFUN <- function(x){sprintf("%.1f", x)}
#'   
#'   # gms
#'   # PROPORTION GERMINATED
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_gMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_germ_1 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS+poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_germ_2 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~poll_treat_factor+(1|session),family=binomial)
#'   mod_gMS_germ_3 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS+(1|session),family=binomial)
#'   anova(mod_gMS_germ_0,mod_gMS_germ_1) # inter
#'   anova(mod_gMS_germ_1,mod_gMS_germ_2) # MS
#'   anova(mod_gMS_germ_1,mod_gMS_germ_3) # ttt
#'   gMS_germ_low_est <- summary(mod_gMS_germ_0)$coefficients["gMS","Estimate"] # low
#'   gMS_germ_low_se <- summary(mod_gMS_germ_0)$coefficients["gMS","Std. Error"] # low
#'   gMS_germ_low_pval <- summary(mod_gMS_germ_0)$coefficients["gMS","Pr(>|z|)"] # low
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_gMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
#'   gMS_germ_med_est <- summary(mod_gMS_germ_0)$coefficients["gMS","Estimate"] # med
#'   gMS_germ_med_se <- summary(mod_gMS_germ_0)$coefficients["gMS","Std. Error"]# med
#'   gMS_germ_med_pval <- summary(mod_gMS_germ_0)$coefficients["gMS","Pr(>|z|)"] # med
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_gMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
#'   gMS_germ_high_est <- summary(mod_gMS_germ_0)$coefficients["gMS","Estimate"] # high
#'   gMS_germ_high_se <- summary(mod_gMS_germ_0)$coefficients["gMS","Std. Error"] # high
#'   gMS_germ_high_pval <- summary(mod_gMS_germ_0)$coefficients["gMS","Pr(>|z|)"]# high
#'   
#'   tab_gMS_germ <- tibble(ttt=c("low","medium","high"),
#'                          est=c(gMS_germ_low_est,gMS_germ_med_est,gMS_germ_high_est),
#'                          se=c(gMS_germ_low_se,gMS_germ_med_se,gMS_germ_high_se),
#'                          pval=c(gMS_germ_low_pval,gMS_germ_med_pval,gMS_germ_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_gMS_germ_fem <- ggplot(tab_gMS_germ, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     xlab("") +
#'     scale_y_continuous(limits=c(-0.5,2.3),labels=scaleFUN) +
#'     ylab("\nGermination rate in relation to\ngenetic number of mates")  + 
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_gMS_germ_fem
#'   ggsave("figures/plot_gMS_germ_fem.jpeg",plot_gMS_germ_fem,width=4,height=4, device = png,bg="white")
#'   
#'   # oms
#'   
#'   # PROPORTION GERMINATED
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_oMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_germ_1 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS+poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_germ_2 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~poll_treat_factor+(1|session),family=binomial)
#'   mod_oMS_germ_3 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS+(1|session),family=binomial)
#'   
#'   anova(mod_oMS_germ_0,mod_oMS_germ_1) # inter
#'   anova(mod_oMS_germ_1,mod_oMS_germ_2) # MS
#'   anova(mod_oMS_germ_1,mod_oMS_germ_3) # ttt
#'   oMS_germ_low_est <- summary(mod_oMS_germ_0)$coefficients["oMS","Estimate"] # low
#'   oMS_germ_low_se <- summary(mod_oMS_germ_0)$coefficients["oMS","Std. Error"] # low
#'   oMS_germ_low_pval <- summary(mod_oMS_germ_0)$coefficients["oMS","Pr(>|z|)"] # low
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_oMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
#'   oMS_germ_med_est <- summary(mod_oMS_germ_0)$coefficients["oMS","Estimate"] # med
#'   oMS_germ_med_se <- summary(mod_oMS_germ_0)$coefficients["oMS","Std. Error"]# med
#'   oMS_germ_med_pval <- summary(mod_oMS_germ_0)$coefficients["oMS","Pr(>|z|)"] # med
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_oMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
#'   oMS_germ_high_est <- summary(mod_oMS_germ_0)$coefficients["oMS","Estimate"] # high
#'   oMS_germ_high_se <- summary(mod_oMS_germ_0)$coefficients["oMS","Std. Error"] # high
#'   oMS_germ_high_pval <- summary(mod_oMS_germ_0)$coefficients["oMS","Pr(>|z|)"]# high
#'   
#'   tab_oMS_germ <- tibble(ttt=c("low","medium","high"),
#'                          est=c(oMS_germ_low_est,oMS_germ_med_est,oMS_germ_high_est),
#'                          se=c(oMS_germ_low_se,oMS_germ_med_se,oMS_germ_high_se),
#'                          pval=c(oMS_germ_low_pval,oMS_germ_med_pval,oMS_germ_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_oMS_germ_fem <- ggplot(tab_oMS_germ, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     xlab("") +
#'     scale_y_continuous(limits=c(-0.7,0.75),labels=scaleFUN) +
#'     ylab("\nGermination rate in relation to\nobservational number of mates")  + 
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_oMS_germ_fem
#'   ggsave("figures/plot_oMS_germ_fem.jpeg",plot_oMS_germ_fem,width=4,height=4, device = png,bg="white")
#'   
#'   
#'   return(list(plot_gms_germ = plot_gMS_germ_fem,
#'               plot_oms_germ = plot_oMS_germ_fem))
#' }
#' 
#' #' Get weight
#' #'
#' #' @description 
#' #'
#' #' @param 
#' #'
#' #' @return 
#' #' 
#' #' @import dplyr
#' #' 
#' #' @export
#' 
#' 
#' get_weight <- function(data_id_fem, common_theme){
#'   
#'   scaleFUN <- function(x){sprintf("%.1f", x)}
#'   
#'   # gms
#'   # MEAN SEED WEIGHT
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_gMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS*poll_treat_factor+(1|session))
#'   mod_gMS_weight_1 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS+poll_treat_factor+(1|session))
#'   mod_gMS_weight_2 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~poll_treat_factor+(1|session))
#'   mod_gMS_weight_3 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS+(1|session))
#'   anova(mod_gMS_weight_0,mod_gMS_weight_1) # inter
#'   anova(mod_gMS_weight_1,mod_gMS_weight_2) # MS
#'   anova(mod_gMS_weight_1,mod_gMS_weight_3) # ttt
#'   gMS_weight_low_est <- summary(mod_gMS_weight_0)$coefficients["gMS","Estimate"] # low
#'   gMS_weight_low_se <- summary(mod_gMS_weight_0)$coefficients["gMS","Std. Error"] # low
#'   gMS_weight_low_pval <- summary(mod_gMS_weight_0)$coefficients["gMS","t value"] # low
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_gMS_weight_0 <- lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS*poll_treat_factor+(1|session))
#'   gMS_weight_med_est <- summary(mod_gMS_weight_0)$coefficients["gMS","Estimate"] # med
#'   gMS_weight_med_se <- summary(mod_gMS_weight_0)$coefficients["gMS","Std. Error"]# med
#'   gMS_weight_med_pval <- summary(mod_gMS_weight_0)$coefficients["gMS","t value"] # med
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_gMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS*poll_treat_factor+(1|session))
#'   gMS_weight_high_est <- summary(mod_gMS_weight_0)$coefficients["gMS","Estimate"] # high
#'   gMS_weight_high_se <- summary(mod_gMS_weight_0)$coefficients["gMS","Std. Error"] # high
#'   gMS_weight_high_pval <- summary(mod_gMS_weight_0)$coefficients["gMS","t value"]# high
#'   
#'   tab_gMS_weight <- tibble(ttt=c("low","medium","high"),
#'                            est=c(gMS_weight_low_est,gMS_weight_med_est,gMS_weight_high_est),
#'                            se=c(gMS_weight_low_se,gMS_weight_med_se,gMS_weight_high_se),
#'                            pval=c(gMS_weight_low_pval,gMS_weight_med_pval,gMS_weight_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_gMS_weight_fem <- ggplot(tab_gMS_weight, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     xlab("") +
#'     ylab("\nSeed weight in relation to\ngenetic number of mates")  + 
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_gMS_weight_fem
#'   ggsave("figures/plot_gMS_weight_fem.jpeg",plot_gMS_weight_fem,width=4,height=4, device = png,bg="white")
#'   
#'   # oms
#'   
#'   # MEAN SEED WEIGHT
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
#'   # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
#'   mod_oMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS*poll_treat_factor+(1|session))
#'   mod_oMS_weight_1 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS+poll_treat_factor+(1|session))
#'   mod_oMS_weight_2 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~poll_treat_factor+(1|session))
#'   mod_oMS_weight_3 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS+(1|session))
#'   anova(mod_oMS_weight_0,mod_oMS_weight_1) # inter
#'   anova(mod_oMS_weight_1,mod_oMS_weight_2) # MS
#'   anova(mod_oMS_weight_1,mod_oMS_weight_3) # ttt
#'   oMS_weight_low_est <- summary(mod_oMS_weight_0)$coefficients["oMS","Estimate"] # low
#'   oMS_weight_low_se <- summary(mod_oMS_weight_0)$coefficients["oMS","Std. Error"] # low
#'   oMS_weight_low_pval <- summary(mod_oMS_weight_0)$coefficients["oMS","t value"] # low
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
#'   mod_oMS_weight_0 <- lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS*poll_treat_factor+(1|session))
#'   oMS_weight_med_est <- summary(mod_oMS_weight_0)$coefficients["oMS","Estimate"] # med
#'   oMS_weight_med_se <- summary(mod_oMS_weight_0)$coefficients["oMS","Std. Error"]# med
#'   oMS_weight_med_pval <- summary(mod_oMS_weight_0)$coefficients["oMS","t value"] # med
#'   data_id_fem <- data_id_fem %>%
#'     mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
#'   mod_oMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS*poll_treat_factor+(1|session))
#'   oMS_weight_high_est <- summary(mod_oMS_weight_0)$coefficients["oMS","Estimate"] # high
#'   oMS_weight_high_se <- summary(mod_oMS_weight_0)$coefficients["oMS","Std. Error"] # high
#'   oMS_weight_high_pval <- summary(mod_oMS_weight_0)$coefficients["oMS","t value"]# high
#'   
#'   tab_oMS_weight <- tibble(ttt=c("low","medium","high"),
#'                            est=c(oMS_weight_low_est,oMS_weight_med_est,oMS_weight_high_est),
#'                            se=c(oMS_weight_low_se,oMS_weight_med_se,oMS_weight_high_se),
#'                            pval=c(oMS_weight_low_pval,oMS_weight_med_pval,oMS_weight_high_pval)) %>%
#'     mutate(est_low=est-1.96*se,
#'            est_up=est+1.96*se) %>%
#'     mutate(pval_stars=case_when(pval<0.001~"***",
#'                                 pval<0.01~"**",
#'                                 pval<0.05~"*",
#'                                 TRUE~"")) %>%
#'     mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
#'   
#'   plot_oMS_weight_fem <- ggplot(tab_oMS_weight, aes(x=ttt, y=est, fill=ttt)) +
#'     geom_hline(yintercept=0,color="gray60",linetype="dashed")+
#'     geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
#'     geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
#'     scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
#'     ggthemes::theme_calc(base_family = "sans") +
#'     common_theme +
#'     xlab("") +
#'     ylab("\nSeed weight in relation to\nobservational number of mates")  + 
#'     scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
#'     geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
#'   plot_oMS_weight_fem
#'   ggsave("figures/plot_oMS_weight_fem.jpeg",plot_oMS_weight_fem,width=4,height=4, device = png,bg="white")
#'   
#'   return(list(plot_gms_weight = plot_gMS_weight_fem,
#'               plot_oms_weight = plot_oMS_weight_fem))
#' }