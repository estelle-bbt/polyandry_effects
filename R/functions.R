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
    mutate(poll_treat_factor=as.factor(case_when(poll_treat==1~"low",
                                                 poll_treat==2~"medium",
                                                 TRUE~"high"))) |>
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high")))
  
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

load_data_both_sexes <- function(file_path, c = "10"){
  
  data_both_sexes <- read.table(file_path,head=T) |>
    mutate(type=as.factor(type),
           poll_treat_factor=as.factor(poll_treat_factor)) |>
    mutate(type=relevel(type, ref = "fem"),
           poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high"))) |>
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

load_data_id <- function(file_path){
  
  data_id <- read.table(file_path,head=T) %>%
    mutate(poll_treat_factor=as.factor(case_when(poll_treat==1~"low",
                                                 poll_treat==2~"medium",
                                                 TRUE~"high"))) %>%
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high"))) %>%
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
  # list_filt <- data_id %>%
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
    left_join(data_id %>% select(ID_full, prop_self)) 
  
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
    left_join(data_id %>% select(ID_full, prop_self))
  
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
    left_join(data_id %>% select(ID_full, prop_self))
  
  
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

#' Get common theme for ggplot
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
get_common_theme <- function(){
  
  common_theme <- theme(axis.text.x = element_text(size=10,angle=45,hjust=1),
                        axis.text.y = element_text(size=8,margin = margin(t = 0, r = 5, b = 0, l = 5)),
                        legend.position = "none",
                        plot.background = element_blank(),
                        panel.grid.major=element_blank(),
                        plot.margin = margin(t=5, r=5, b=5, l=5),
                        strip.background = element_rect(colour="black", 
                                                        linewidth=1.5, linetype="solid"),
                        panel.background = element_rect(colour="black", fill="white",linewidth=1),
                        axis.title.y = element_text(size=12))
  return(common_theme)
}

#' Get proportion of non-aborted
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


get_prop_non_abort <- function(data_proxy_av, common_theme){
  
  scaleFUN <- function(x){sprintf("%.1f", x)}
  
  # gms
  
  data_proxy_av <- data_proxy_av %>%
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_gMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS*poll_treat_factor+(1|session),family=binomial)
  mod_gMS_av_1 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS+poll_treat_factor+(1|session),family=binomial)
  mod_gMS_av_2 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~poll_treat_factor+(1|session),family=binomial)
  mod_gMS_av_3 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS+(1|session),family=binomial)
  anova(mod_gMS_av_0,mod_gMS_av_1) # inter
  anova(mod_gMS_av_1,mod_gMS_av_2) # MS
  anova(mod_gMS_av_1,mod_gMS_av_3) # ttt
  gMS_av_low_est <- summary(mod_gMS_av_0)$coefficients["gMS","Estimate"] # low
  gMS_av_low_se <- summary(mod_gMS_av_0)$coefficients["gMS","Std. Error"] # low
  gMS_av_low_pval <- summary(mod_gMS_av_0)$coefficients["gMS","Pr(>|z|)"] # low
  data_proxy_av <- data_proxy_av %>%
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_gMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS*poll_treat_factor+(1|session),family=binomial)
  gMS_av_med_est <- summary(mod_gMS_av_0)$coefficients["gMS","Estimate"] # med
  gMS_av_med_se <- summary(mod_gMS_av_0)$coefficients["gMS","Std. Error"]# med
  gMS_av_med_pval <- summary(mod_gMS_av_0)$coefficients["gMS","Pr(>|z|)"] # med
  data_proxy_av <- data_proxy_av %>%
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_gMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(gMS))),cbind(nbGr_sum,nbAv_sum)~gMS*poll_treat_factor+(1|session),family=binomial)
  gMS_av_high_est <- summary(mod_gMS_av_0)$coefficients["gMS","Estimate"] # high
  gMS_av_high_se <- summary(mod_gMS_av_0)$coefficients["gMS","Std. Error"] # high
  gMS_av_high_pval <- summary(mod_gMS_av_0)$coefficients["gMS","Pr(>|z|)"]# high
  
  tab_gMS_av <- tibble(ttt=c("low","medium","high"),
                       est=c(gMS_av_low_est,gMS_av_med_est,gMS_av_high_est),
                       se=c(gMS_av_low_se,gMS_av_med_se,gMS_av_high_se),
                       pval=c(gMS_av_low_pval,gMS_av_med_pval,gMS_av_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=forcats::fct_relevel(ttt, c("low","medium","high")))
  
  plot_gMS_av_fem <- ggplot(tab_gMS_av, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    xlab("") +
    ylab("\nProportion of non-aborted seeds \nin relation to genetic number of mates")  + 
    ylim(-0.3,1.1) +
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_gMS_av_fem
  ggsave("figures/plot_gMS_av_fem.jpeg",plot_gMS_av_fem,width=4,height=4, device = png,bg="white")
  
  # oms
  
  # PROPORTION NON-ABORTED
  data_proxy_av <- data_proxy_av %>%
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_oMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS*poll_treat_factor+(1|session),family=binomial)
  mod_oMS_av_1 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS+poll_treat_factor+(1|session),family=binomial)
  mod_oMS_av_2 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~poll_treat_factor+(1|session),family=binomial)
  mod_oMS_av_3 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS+(1|session),family=binomial)
  
  anova(mod_oMS_av_0,mod_oMS_av_1) # inter
  anova(mod_oMS_av_1,mod_oMS_av_2) # MS
  anova(mod_oMS_av_1,mod_oMS_av_3) # ttt
  oMS_av_low_est <- summary(mod_oMS_av_0)$coefficients["oMS","Estimate"] # low
  oMS_av_low_se <- summary(mod_oMS_av_0)$coefficients["oMS","Std. Error"] # low
  oMS_av_low_pval <- summary(mod_oMS_av_0)$coefficients["oMS","Pr(>|z|)"] # low
  data_proxy_av <- data_proxy_av %>%
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_oMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS*poll_treat_factor+(1|session),family=binomial)
  oMS_av_med_est <- summary(mod_oMS_av_0)$coefficients["oMS","Estimate"] # med
  oMS_av_med_se <- summary(mod_oMS_av_0)$coefficients["oMS","Std. Error"]# med
  oMS_av_med_pval <- summary(mod_oMS_av_0)$coefficients["oMS","Pr(>|z|)"] # med
  data_proxy_av <- data_proxy_av %>%
    mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_oMS_av_0 <-  lme4::glmer(data=(data_proxy_av %>% filter(!is.na(oMS))),cbind(nbGr_sum,nbAv_sum)~oMS*poll_treat_factor+(1|session),family=binomial)
  oMS_av_high_est <- summary(mod_oMS_av_0)$coefficients["oMS","Estimate"] # high
  oMS_av_high_se <- summary(mod_oMS_av_0)$coefficients["oMS","Std. Error"] # high
  oMS_av_high_pval <- summary(mod_oMS_av_0)$coefficients["oMS","Pr(>|z|)"]# high
  
  tab_oMS_av <- tibble(ttt=c("low","medium","high"),
                       est=c(oMS_av_low_est,oMS_av_med_est,oMS_av_high_est),
                       se=c(oMS_av_low_se,oMS_av_med_se,oMS_av_high_se),
                       pval=c(oMS_av_low_pval,oMS_av_med_pval,oMS_av_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=forcats::fct_relevel(ttt, c("low","medium","high")))
  
  plot_oMS_av_fem <- ggplot(tab_oMS_av, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    xlab("") +
    scale_y_continuous(limits=c(-0.6,0.8),labels=scaleFUN) +
    ylab("\nProportion of non-aborted seeds \nin relation to observational number of mates")  + 
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_oMS_av_fem
  ggsave("figures/plot_oMS_av_fem.jpeg",plot_oMS_av_fem,width=4,height=4, device = png,bg="white")
  
  
  
  return(list(plot_gms_av = plot_gMS_av_fem,
              plot_oms_av = plot_oMS_av_fem))
}

#' Get seed-set
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


get_seed_set <- function(data_proxy_ov, common_theme){
  
  scaleFUN <- function(x){sprintf("%.1f", x)}
  
  # gms
  
  data_proxy_ov <- data_proxy_ov %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_gMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
  mod_gMS_ss_1 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS+poll_treat_factor+(1|session),family=binomial)
  mod_gMS_ss_2 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~poll_treat_factor+(1|session),family=binomial)
  mod_gMS_ss_3 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS+(1|session),family=binomial)
  anova(mod_gMS_ss_0,mod_gMS_ss_1) # inter
  anova(mod_gMS_ss_1,mod_gMS_ss_2) # MS
  anova(mod_gMS_ss_1,mod_gMS_ss_3) # ttt
  gMS_ss_low_est <- summary(mod_gMS_ss_0)$coefficients["gMS","Estimate"] # low
  gMS_ss_low_se <- summary(mod_gMS_ss_0)$coefficients["gMS","Std. Error"] # low
  gMS_ss_low_pval <- summary(mod_gMS_ss_0)$coefficients["gMS","Pr(>|z|)"] # low
  data_proxy_ov <- data_proxy_ov %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_gMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
  gMS_ss_med_est <- summary(mod_gMS_ss_0)$coefficients["gMS","Estimate"] # med
  gMS_ss_med_se <- summary(mod_gMS_ss_0)$coefficients["gMS","Std. Error"]# med
  gMS_ss_med_pval <- summary(mod_gMS_ss_0)$coefficients["gMS","Pr(>|z|)"] # med
  data_proxy_ov <- data_proxy_ov %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_gMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(gMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
  gMS_ss_high_est <- summary(mod_gMS_ss_0)$coefficients["gMS","Estimate"] # high
  gMS_ss_high_se <- summary(mod_gMS_ss_0)$coefficients["gMS","Std. Error"] # high
  gMS_ss_high_pval <- summary(mod_gMS_ss_0)$coefficients["gMS","Pr(>|z|)"]# high
  
  tab_gMS_ss <- tibble(ttt=c("low","medium","high"),
                       est=c(gMS_ss_low_est,gMS_ss_med_est,gMS_ss_high_est),
                       se=c(gMS_ss_low_se,gMS_ss_med_se,gMS_ss_high_se),
                       pval=c(gMS_ss_low_pval,gMS_ss_med_pval,gMS_ss_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
  
  plot_gMS_ss_fem <- ggplot(tab_gMS_ss, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    ylim(0,0.5) +
    xlab("") +
    ylab("\nSeed-set in relation to\ngenetic number of mates")  + 
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_gMS_ss_fem
  ggsave("figures/plot_gMS_ss_fem.jpeg",plot_gMS_ss_fem,width=4,height=4, device = png,bg="white")
  
  # oms
  
  data_proxy_ov <- data_proxy_ov %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_oMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
  mod_oMS_ss_1 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS+poll_treat_factor+(1|session),family=binomial)
  mod_oMS_ss_2 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~poll_treat_factor+(1|session),family=binomial)
  mod_oMS_ss_3 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS+(1|session),family=binomial)
  
  anova(mod_oMS_ss_0,mod_oMS_ss_1) # inter
  anova(mod_oMS_ss_1,mod_oMS_ss_2) # MS
  anova(mod_oMS_ss_1,mod_oMS_ss_3) # ttt
  oMS_ss_low_est <- summary(mod_oMS_ss_0)$coefficients["oMS","Estimate"] # low
  oMS_ss_low_se <- summary(mod_oMS_ss_0)$coefficients["oMS","Std. Error"] # low
  oMS_ss_low_pval <- summary(mod_oMS_ss_0)$coefficients["oMS","Pr(>|z|)"] # low
  data_proxy_ov <- data_proxy_ov %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_oMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
  oMS_ss_med_est <- summary(mod_oMS_ss_0)$coefficients["oMS","Estimate"] # med
  oMS_ss_med_se <- summary(mod_oMS_ss_0)$coefficients["oMS","Std. Error"]# med
  oMS_ss_med_pval <- summary(mod_oMS_ss_0)$coefficients["oMS","Pr(>|z|)"] # med
  data_proxy_ov <- data_proxy_ov %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_oMS_ss_0 <-  lme4::glmer(data=(data_proxy_ov %>% filter(!is.na(oMS))),cbind(nbGr_sum,(nbOv_sum-nbGr_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
  oMS_ss_high_est <- summary(mod_oMS_ss_0)$coefficients["oMS","Estimate"] # high
  oMS_ss_high_se <- summary(mod_oMS_ss_0)$coefficients["oMS","Std. Error"] # high
  oMS_ss_high_pval <- summary(mod_oMS_ss_0)$coefficients["oMS","Pr(>|z|)"]# high
  
  tab_oMS_ss <- tibble(ttt=c("low","medium","high"),
                       est=c(oMS_ss_low_est,oMS_ss_med_est,oMS_ss_high_est),
                       se=c(oMS_ss_low_se,oMS_ss_med_se,oMS_ss_high_se),
                       pval=c(oMS_ss_low_pval,oMS_ss_med_pval,oMS_ss_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
  
  plot_oMS_ss_fem <- ggplot(tab_oMS_ss, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    xlab("") +
    scale_y_continuous(limits=c(-0.15,0.4),labels=scaleFUN) +
    ylab("\nSeed-set in relation to\nobservational number of mates")  + 
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_oMS_ss_fem
  ggsave("figures/plot_oMS_ss_fem.jpeg",plot_oMS_ss_fem,width=4,height=4, device = png,bg="white")
  
  
  return(list(plot_gms_ss = plot_gMS_ss_fem,
              plot_oms_ss = plot_oMS_ss_fem))
}

#' Get proportion germinated seeds
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


get_prop_germ <- function(data_id_fem, common_theme){
  
  scaleFUN <- function(x){sprintf("%.1f", x)}
  
  # gms
  # PROPORTION GERMINATED
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_gMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
  mod_gMS_germ_1 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS+poll_treat_factor+(1|session),family=binomial)
  mod_gMS_germ_2 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~poll_treat_factor+(1|session),family=binomial)
  mod_gMS_germ_3 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS+(1|session),family=binomial)
  anova(mod_gMS_germ_0,mod_gMS_germ_1) # inter
  anova(mod_gMS_germ_1,mod_gMS_germ_2) # MS
  anova(mod_gMS_germ_1,mod_gMS_germ_3) # ttt
  gMS_germ_low_est <- summary(mod_gMS_germ_0)$coefficients["gMS","Estimate"] # low
  gMS_germ_low_se <- summary(mod_gMS_germ_0)$coefficients["gMS","Std. Error"] # low
  gMS_germ_low_pval <- summary(mod_gMS_germ_0)$coefficients["gMS","Pr(>|z|)"] # low
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_gMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
  gMS_germ_med_est <- summary(mod_gMS_germ_0)$coefficients["gMS","Estimate"] # med
  gMS_germ_med_se <- summary(mod_gMS_germ_0)$coefficients["gMS","Std. Error"]# med
  gMS_germ_med_pval <- summary(mod_gMS_germ_0)$coefficients["gMS","Pr(>|z|)"] # med
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_gMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(gMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~gMS*poll_treat_factor+(1|session),family=binomial)
  gMS_germ_high_est <- summary(mod_gMS_germ_0)$coefficients["gMS","Estimate"] # high
  gMS_germ_high_se <- summary(mod_gMS_germ_0)$coefficients["gMS","Std. Error"] # high
  gMS_germ_high_pval <- summary(mod_gMS_germ_0)$coefficients["gMS","Pr(>|z|)"]# high
  
  tab_gMS_germ <- tibble(ttt=c("low","medium","high"),
                         est=c(gMS_germ_low_est,gMS_germ_med_est,gMS_germ_high_est),
                         se=c(gMS_germ_low_se,gMS_germ_med_se,gMS_germ_high_se),
                         pval=c(gMS_germ_low_pval,gMS_germ_med_pval,gMS_germ_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
  
  plot_gMS_germ_fem <- ggplot(tab_gMS_germ, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    xlab("") +
    scale_y_continuous(limits=c(-0.5,2.3),labels=scaleFUN) +
    ylab("\nGermination rate in relation to\ngenetic number of mates")  + 
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_gMS_germ_fem
  ggsave("figures/plot_gMS_germ_fem.jpeg",plot_gMS_germ_fem,width=4,height=4, device = png,bg="white")
  
  # oms
  
  # PROPORTION GERMINATED
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_oMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
  mod_oMS_germ_1 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS+poll_treat_factor+(1|session),family=binomial)
  mod_oMS_germ_2 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~poll_treat_factor+(1|session),family=binomial)
  mod_oMS_germ_3 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS+(1|session),family=binomial)
  
  anova(mod_oMS_germ_0,mod_oMS_germ_1) # inter
  anova(mod_oMS_germ_1,mod_oMS_germ_2) # MS
  anova(mod_oMS_germ_1,mod_oMS_germ_3) # ttt
  oMS_germ_low_est <- summary(mod_oMS_germ_0)$coefficients["oMS","Estimate"] # low
  oMS_germ_low_se <- summary(mod_oMS_germ_0)$coefficients["oMS","Std. Error"] # low
  oMS_germ_low_pval <- summary(mod_oMS_germ_0)$coefficients["oMS","Pr(>|z|)"] # low
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_oMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
  oMS_germ_med_est <- summary(mod_oMS_germ_0)$coefficients["oMS","Estimate"] # med
  oMS_germ_med_se <- summary(mod_oMS_germ_0)$coefficients["oMS","Std. Error"]# med
  oMS_germ_med_pval <- summary(mod_oMS_germ_0)$coefficients["oMS","Pr(>|z|)"] # med
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_oMS_germ_0 <-  lme4::glmer(data=(data_id_fem %>% filter(!is.na(oMS))),cbind(seed_germ_sum,(seed_semis_sum-seed_germ_sum))~oMS*poll_treat_factor+(1|session),family=binomial)
  oMS_germ_high_est <- summary(mod_oMS_germ_0)$coefficients["oMS","Estimate"] # high
  oMS_germ_high_se <- summary(mod_oMS_germ_0)$coefficients["oMS","Std. Error"] # high
  oMS_germ_high_pval <- summary(mod_oMS_germ_0)$coefficients["oMS","Pr(>|z|)"]# high
  
  tab_oMS_germ <- tibble(ttt=c("low","medium","high"),
                         est=c(oMS_germ_low_est,oMS_germ_med_est,oMS_germ_high_est),
                         se=c(oMS_germ_low_se,oMS_germ_med_se,oMS_germ_high_se),
                         pval=c(oMS_germ_low_pval,oMS_germ_med_pval,oMS_germ_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
  
  plot_oMS_germ_fem <- ggplot(tab_oMS_germ, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    xlab("") +
    scale_y_continuous(limits=c(-0.7,0.75),labels=scaleFUN) +
    ylab("\nGermination rate in relation to\nobservational number of mates")  + 
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_oMS_germ_fem
  ggsave("figures/plot_oMS_germ_fem.jpeg",plot_oMS_germ_fem,width=4,height=4, device = png,bg="white")
  
  
  return(list(plot_gms_germ = plot_gMS_germ_fem,
              plot_oms_germ = plot_oMS_germ_fem))
}

#' Get weight
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


get_weight <- function(data_id_fem, common_theme){
  
  scaleFUN <- function(x){sprintf("%.1f", x)}
  
  # gms
  # MEAN SEED WEIGHT
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_gMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS*poll_treat_factor+(1|session))
  mod_gMS_weight_1 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS+poll_treat_factor+(1|session))
  mod_gMS_weight_2 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~poll_treat_factor+(1|session))
  mod_gMS_weight_3 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS+(1|session))
  anova(mod_gMS_weight_0,mod_gMS_weight_1) # inter
  anova(mod_gMS_weight_1,mod_gMS_weight_2) # MS
  anova(mod_gMS_weight_1,mod_gMS_weight_3) # ttt
  gMS_weight_low_est <- summary(mod_gMS_weight_0)$coefficients["gMS","Estimate"] # low
  gMS_weight_low_se <- summary(mod_gMS_weight_0)$coefficients["gMS","Std. Error"] # low
  gMS_weight_low_pval <- summary(mod_gMS_weight_0)$coefficients["gMS","t value"] # low
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_gMS_weight_0 <- lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS*poll_treat_factor+(1|session))
  gMS_weight_med_est <- summary(mod_gMS_weight_0)$coefficients["gMS","Estimate"] # med
  gMS_weight_med_se <- summary(mod_gMS_weight_0)$coefficients["gMS","Std. Error"]# med
  gMS_weight_med_pval <- summary(mod_gMS_weight_0)$coefficients["gMS","t value"] # med
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_gMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(gMS))),poids_sd_mean~gMS*poll_treat_factor+(1|session))
  gMS_weight_high_est <- summary(mod_gMS_weight_0)$coefficients["gMS","Estimate"] # high
  gMS_weight_high_se <- summary(mod_gMS_weight_0)$coefficients["gMS","Std. Error"] # high
  gMS_weight_high_pval <- summary(mod_gMS_weight_0)$coefficients["gMS","t value"]# high
  
  tab_gMS_weight <- tibble(ttt=c("low","medium","high"),
                           est=c(gMS_weight_low_est,gMS_weight_med_est,gMS_weight_high_est),
                           se=c(gMS_weight_low_se,gMS_weight_med_se,gMS_weight_high_se),
                           pval=c(gMS_weight_low_pval,gMS_weight_med_pval,gMS_weight_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
  
  plot_gMS_weight_fem <- ggplot(tab_gMS_weight, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    xlab("") +
    ylab("\nSeed weight in relation to\ngenetic number of mates")  + 
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_gMS_weight_fem
  ggsave("figures/plot_gMS_weight_fem.jpeg",plot_gMS_weight_fem,width=4,height=4, device = png,bg="white")
  
  # oms
  
  # MEAN SEED WEIGHT
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("low","medium","high")))
  # to apply likelihood ratio test, models must applied to the same data set -> always remove row with NA
  mod_oMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS*poll_treat_factor+(1|session))
  mod_oMS_weight_1 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS+poll_treat_factor+(1|session))
  mod_oMS_weight_2 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~poll_treat_factor+(1|session))
  mod_oMS_weight_3 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS+(1|session))
  anova(mod_oMS_weight_0,mod_oMS_weight_1) # inter
  anova(mod_oMS_weight_1,mod_oMS_weight_2) # MS
  anova(mod_oMS_weight_1,mod_oMS_weight_3) # ttt
  oMS_weight_low_est <- summary(mod_oMS_weight_0)$coefficients["oMS","Estimate"] # low
  oMS_weight_low_se <- summary(mod_oMS_weight_0)$coefficients["oMS","Std. Error"] # low
  oMS_weight_low_pval <- summary(mod_oMS_weight_0)$coefficients["oMS","t value"] # low
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("medium","high","low")))
  mod_oMS_weight_0 <- lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS*poll_treat_factor+(1|session))
  oMS_weight_med_est <- summary(mod_oMS_weight_0)$coefficients["oMS","Estimate"] # med
  oMS_weight_med_se <- summary(mod_oMS_weight_0)$coefficients["oMS","Std. Error"]# med
  oMS_weight_med_pval <- summary(mod_oMS_weight_0)$coefficients["oMS","t value"] # med
  data_id_fem <- data_id_fem %>%
    mutate(poll_treat_factor=fct_relevel(poll_treat_factor, c("high","low","medium")))
  mod_oMS_weight_0 <-  lme4::lmer(data=(data_id_fem %>% filter(!is.na(oMS))),poids_sd_mean~oMS*poll_treat_factor+(1|session))
  oMS_weight_high_est <- summary(mod_oMS_weight_0)$coefficients["oMS","Estimate"] # high
  oMS_weight_high_se <- summary(mod_oMS_weight_0)$coefficients["oMS","Std. Error"] # high
  oMS_weight_high_pval <- summary(mod_oMS_weight_0)$coefficients["oMS","t value"]# high
  
  tab_oMS_weight <- tibble(ttt=c("low","medium","high"),
                           est=c(oMS_weight_low_est,oMS_weight_med_est,oMS_weight_high_est),
                           se=c(oMS_weight_low_se,oMS_weight_med_se,oMS_weight_high_se),
                           pval=c(oMS_weight_low_pval,oMS_weight_med_pval,oMS_weight_high_pval)) %>%
    mutate(est_low=est-1.96*se,
           est_up=est+1.96*se) %>%
    mutate(pval_stars=case_when(pval<0.001~"***",
                                pval<0.01~"**",
                                pval<0.05~"*",
                                TRUE~"")) %>%
    mutate(ttt=fct_relevel(ttt, c("low","medium","high")))
  
  plot_oMS_weight_fem <- ggplot(tab_oMS_weight, aes(x=ttt, y=est, fill=ttt)) +
    geom_hline(yintercept=0,color="gray60",linetype="dashed")+
    geom_errorbar(aes(ymin=est_low, ymax=est_up,color=ttt), stat="identity",width=0,linewidth=3,lineend="butt") +
    geom_point(size=4,shape=21,color="gray30",fill="white",stroke=1) +
    scale_fill_manual(values=c("#DCBC93","#C89656","#81675F")) +
    scale_color_manual(values=c("#DCBC93","#C89656","#81675F")) +
    ggthemes::theme_calc(base_family = "sans") +
    common_theme +
    xlab("") +
    ylab("\nSeed weight in relation to\nobservational number of mates")  + 
    scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
    geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
  plot_oMS_weight_fem
  ggsave("figures/plot_oMS_weight_fem.jpeg",plot_oMS_weight_fem,width=4,height=4, device = png,bg="white")
  
  return(list(plot_gms_weight = plot_gMS_weight_fem,
              plot_oms_weight = plot_oMS_weight_fem))
}