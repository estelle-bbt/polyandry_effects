library(tidyverse)
targets::tar_load(data_flower)
targets::tar_load(data_both_sexes)

dt <- data_flower %>%
  left_join(data_both_sexes %>%
              filter(type=="fem") %>%
              select(ID_full,gMS,nb_part_ID_out_co10), by = join_by(ID_full)) %>% # to keep the same dataset
  dplyr::mutate( # if convergence problem, center-scale continuous
    nb_part_ID_out_co10_c = scale(nb_part_ID_out_co10, center = TRUE, scale = TRUE),
    import_nb_part_ID_out_co10_c = scale(import_nb_part_ID_out_co10, center = TRUE, scale = TRUE),
    pl_call = scale(pl, center = TRUE, scale = TRUE)
  )

dt_sum <- dt %>%
  mutate(prop_non_av = nbGr/(nbGr+nbAv),
         prop_germ = nbGr_germ/nbGr_semis,
         seed_set = nbGr_SR/nbOv) %>%
  group_by(poll_treat_factor) %>%
  summarise(mean_prop_non_av = mean(prop_non_av,na.rm=T),
            mean_prop_germ = mean(prop_germ,na.rm=T),
            mean_seed_set = mean(seed_set,na.rm=T),
            mean_weight_seed = mean(poids_sd,na.rm=T))
  
scaleFUN <- function(x){sprintf("%.1f", x)}

# 1) controlling for oms at the flower scale (complete dataset) --------------

# prop avorted
coms_av_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full), family = binomial)
coms_av_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr, nbAv) ~ gMS + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_av_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr, nbAv) ~ poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full), family = binomial)
coms_av_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr, nbAv) ~ gMS + import_nb_part_ID_out_co10  + (1|session:ID_full), family = binomial)
anova(coms_av_gms0,coms_av_gms1)
anova(coms_av_gms1,coms_av_gms2)
anova(coms_av_gms1,coms_av_gms3)

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high")))
coms_av_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
gMS_av_low_est <- summary(coms_av_gms0)$coefficients["gMS","Estimate"] # low
gMS_av_low_se <- summary(coms_av_gms0)$coefficients["gMS","Std. Error"] # low
gMS_av_low_pval <- summary(coms_av_gms0)$coefficients["gMS","Pr(>|z|)"] # low
dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("medium","high","low")))
coms_av_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ gMS * poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
gMS_av_med_est <- summary(coms_av_gms0)$coefficients["gMS","Estimate"] # med
gMS_av_med_se <- summary(coms_av_gms0)$coefficients["gMS","Std. Error"]# med
gMS_av_med_pval <- summary(coms_av_gms0)$coefficients["gMS","Pr(>|z|)"] # med
dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("high","low","medium")))
coms_av_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
gMS_av_high_est <- summary(coms_av_gms0)$coefficients["gMS","Estimate"] # high
gMS_av_high_se <- summary(coms_av_gms0)$coefficients["gMS","Std. Error"] # high
gMS_av_high_pval <- summary(coms_av_gms0)$coefficients["gMS","Pr(>|z|)"]# high

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
  get_common_theme() +
  xlab("") +
  ylab("\nProportion of non-aborted seeds \nin relation to genetic number of mates")  + 
  # ylim(-0.3,1.1) +
  scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
  geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
plot_gMS_av_fem

coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_av_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_av_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_av_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_av_oms0,coms_av_oms1)
anova(coms_av_oms1,coms_av_oms2)
anova(coms_av_oms1,coms_av_oms3)

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("low","medium","high")))
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
oMS_av_low_est <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Estimate"] # low
oMS_av_low_se <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Std. Error"] # low
oMS_av_low_pval <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Pr(>|z|)"] # low
dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("medium","high","low")))
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
oMS_av_med_est <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Estimate"] # med
oMS_av_med_se <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Std. Error"]# med
oMS_av_med_pval <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Pr(>|z|)"] # med
dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, c("high","low","medium")))
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
oMS_av_high_est <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Estimate"] # high
oMS_av_high_se <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Std. Error"] # high
oMS_av_high_pval <- summary(coms_av_oms0)$coefficients["nb_part_ID_out_co10","Pr(>|z|)"]# high

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
  get_common_theme() +
  xlab("") +
  ylab("\nProportion of non-aborted seeds \nin relation to observational number of mates")  + 
  # ylim(-0.3,1.1) +
  scale_x_discrete(labels=c("low" = "LOW", "medium" = "MEDIUM", "high" = "HIGH")) +
  geom_text(aes(y=est_up,label=pval_stars,color=ttt),vjust=-0.2, size=4) 
plot_oMS_av_fem

# prop germinated
coms_germ_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_germ_gms0,coms_germ_gms1)
anova(coms_germ_gms1,coms_germ_gms2)
anova(coms_germ_gms1,coms_germ_gms3)

coms_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_germ_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
anova(coms_germ_oms0,coms_germ_oms1)
anova(coms_germ_oms1,coms_germ_oms2)
anova(coms_germ_oms1,coms_germ_oms3)

# seed-set
coms_ss_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_ss_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_ss_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_ss_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
anova(coms_ss_gms0,coms_ss_gms1)
anova(coms_ss_gms1,coms_ss_gms2)
anova(coms_ss_gms1,coms_ss_gms3)

coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10* poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10+ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor+import_nb_part_ID_out_co10    + (1|session:ID_full),family = binomial)
coms_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10    + (1|session:ID_full),family = binomial)
anova(coms_ss_oms0,coms_ss_oms1)
anova(coms_ss_oms1,coms_ss_oms2)
anova(coms_ss_oms1,coms_ss_oms3)

# seed weight
coms_weight_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+import_nb_part_ID_out_co10  + (1|session:ID_full))
anova(coms_weight_gms0,coms_weight_gms1)
anova(coms_weight_gms1,coms_weight_gms2)
anova(coms_weight_gms1,coms_weight_gms3)

coms_weight_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full))
coms_weight_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10   + (1|session:ID_full))
anova(coms_weight_oms0,coms_weight_oms1)
anova(coms_weight_oms1,coms_weight_oms2)
anova(coms_weight_oms1,coms_weight_oms3)
# 5 effects


# 2) controlling further for pl (reduced dataset) ----------------------------

# prop avorted
comspl_av_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_av_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_av_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_av_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ gMS+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
anova(comspl_av_gms0,comspl_av_gms1)
anova(comspl_av_gms1,comspl_av_gms2)
anova(comspl_av_gms1,comspl_av_gms3)

comspl_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_av_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_av_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_av_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
anova(comspl_av_oms0,comspl_av_oms1)
anova(comspl_av_oms1,comspl_av_oms2)
anova(comspl_av_oms1,comspl_av_oms3)

# prop germinated
comspl_germ_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_germ_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_germ_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_germ_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
anova(comspl_germ_gms0,comspl_germ_gms1)
anova(comspl_germ_gms1,comspl_germ_gms2)
anova(comspl_germ_gms1,comspl_germ_gms3)

comspl_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_germ_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_germ_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor+import_nb_part_ID_out_co10+pl   + (1|session:ID_full),family = binomial)
comspl_germ_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10+pl   + (1|session:ID_full),family = binomial)
anova(comspl_germ_oms0,comspl_germ_oms1)
anova(comspl_germ_oms1,comspl_germ_oms2)
anova(comspl_germ_oms1,comspl_germ_oms3)

# seed-set
comspl_ss_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10+pl   + (1|session:ID_full),family = binomial)
comspl_ss_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10+pl   + (1|session:ID_full),family = binomial)
comspl_ss_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor+import_nb_part_ID_out_co10+pl   + (1|session:ID_full),family = binomial)
comspl_ss_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS+import_nb_part_ID_out_co10+pl   + (1|session:ID_full),family = binomial)
anova(comspl_ss_gms0,comspl_ss_gms1)
anova(comspl_ss_gms1,comspl_ss_gms2)
anova(comspl_ss_gms1,comspl_ss_gms3)

comspl_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10* poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10+ poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full),family = binomial)
comspl_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor+import_nb_part_ID_out_co10+pl    + (1|session:ID_full),family = binomial)
comspl_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10+pl    + (1|session:ID_full),family = binomial)
anova(comspl_ss_oms0,comspl_ss_oms1)
anova(comspl_ss_oms1,comspl_ss_oms2)
anova(comspl_ss_oms1,comspl_ss_oms3)

# seed weight
comspl_weight_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full))
comspl_weight_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full))
comspl_weight_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full))
comspl_weight_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+import_nb_part_ID_out_co10+pl  + (1|session:ID_full))
anova(comspl_weight_gms0,comspl_weight_gms1)
anova(comspl_weight_gms1,comspl_weight_gms2)
anova(comspl_weight_gms1,comspl_weight_gms3)

comspl_weight_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full))
comspl_weight_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10+pl  + (1|session:ID_full))
comspl_weight_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10+pl   + (1|session:ID_full))
comspl_weight_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10+pl   + (1|session:ID_full))
anova(comspl_weight_oms0,comspl_weight_oms1)
anova(comspl_weight_oms1,comspl_weight_oms2)
anova(comspl_weight_oms1,comspl_weight_oms3)
# not any effect with the reduced dataset with pl


# 3) controlling for nothing -------------------------------------------------

# prop avorted
not_av_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ gMS*poll_treat_factor  + (1|session:ID_full),family = binomial)
not_av_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ gMS+poll_treat_factor  + (1|session:ID_full),family = binomial)
not_av_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ poll_treat_factor  + (1|session:ID_full),family = binomial)
not_av_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr,nbAv) ~ gMS  + (1|session:ID_full),family = binomial)
anova(not_av_gms0,not_av_gms1)
anova(not_av_gms1,not_av_gms2)
anova(not_av_gms1,not_av_gms3)

not_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10*poll_treat_factor  + (1|session:ID_full),family = binomial)
not_av_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10+poll_treat_factor  + (1|session:ID_full),family = binomial)
not_av_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ poll_treat_factor  + (1|session:ID_full),family = binomial)
not_av_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(not_av_oms0,not_av_oms1)
anova(not_av_oms1,not_av_oms2)
anova(not_av_oms1,not_av_oms3)

# prop germinated
not_germ_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS*poll_treat_factor  + (1|session:ID_full),family = binomial)
not_germ_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS+poll_treat_factor  + (1|session:ID_full),family = binomial)
not_germ_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor  + (1|session:ID_full),family = binomial)
not_germ_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS  + (1|session:ID_full),family = binomial)
anova(not_germ_gms0,not_germ_gms1)
anova(not_germ_gms1,not_germ_gms2)
anova(not_germ_gms1,not_germ_gms3)

not_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10*poll_treat_factor  + (1|session:ID_full),family = binomial)
not_germ_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10+poll_treat_factor  + (1|session:ID_full),family = binomial)
not_germ_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor   + (1|session:ID_full),family = binomial)
not_germ_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
anova(not_germ_oms0,not_germ_oms1)
anova(not_germ_oms1,not_germ_oms2)
anova(not_germ_oms1,not_germ_oms3)

# seed-set
not_ss_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS*poll_treat_factor   + (1|session:ID_full),family = binomial)
not_ss_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS+poll_treat_factor   + (1|session:ID_full),family = binomial)
not_ss_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor   + (1|session:ID_full),family = binomial)
not_ss_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS   + (1|session:ID_full),family = binomial)
anova(not_ss_gms0,not_ss_gms1)
anova(not_ss_gms1,not_ss_gms2)
anova(not_ss_gms1,not_ss_gms3)

not_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10* poll_treat_factor  + (1|session:ID_full),family = binomial)
not_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10+ poll_treat_factor  + (1|session:ID_full),family = binomial)
not_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor    + (1|session:ID_full),family = binomial)
not_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10    + (1|session:ID_full),family = binomial)
anova(not_ss_oms0,not_ss_oms1)
anova(not_ss_oms1,not_ss_oms2)
anova(not_ss_oms1,not_ss_oms3)

# seed weight
not_weight_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS*poll_treat_factor  + (1|session:ID_full))
not_weight_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+poll_treat_factor  + (1|session:ID_full))
not_weight_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ poll_treat_factor  + (1|session:ID_full))
not_weight_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS  + (1|session:ID_full))
anova(not_weight_gms0,not_weight_gms1)
anova(not_weight_gms1,not_weight_gms2)
anova(not_weight_gms1,not_weight_gms3)

not_weight_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10*poll_treat_factor  + (1|session:ID_full))
not_weight_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10+poll_treat_factor  + (1|session:ID_full))
not_weight_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ poll_treat_factor   + (1|session:ID_full))
not_weight_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10   + (1|session:ID_full))
anova(not_weight_oms0,not_weight_oms1)
anova(not_weight_oms1,not_weight_oms2)
anova(not_weight_oms1,not_weight_oms3)
# 3 effects


# 4) focusing on the flower scale to complete those at the individual scale --------

test <- dt %>% filter(gMS_fem_out > 1)

# prop avorted

# gms - nGenot as covariable
coms_av_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS_fem_out)), cbind(nbGr, nbAv) ~ gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS_fem_out)), cbind(nbGr, nbAv) ~ gMS_fem_out + poll_treat_factor + nGenot  + (1|session:ID_full),family = binomial)
coms_av_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS_fem_out)), cbind(nbGr, nbAv) ~ poll_treat_factor  + nGenot + (1|session:ID_full), family = binomial)
coms_av_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS_fem_out)), cbind(nbGr, nbAv) ~ gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
anova(coms_av_gms0,coms_av_gms1)
anova(coms_av_gms1,coms_av_gms2)
anova(coms_av_gms1,coms_av_gms3)

# oms - pl as covariable
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ import_nb_part_ID_out_co10*poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_av_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ import_nb_part_ID_out_co10+poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_av_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_av_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)), cbind(nbGr,nbAv) ~ import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
anova(coms_av_oms0,coms_av_oms1)
anova(coms_av_oms1,coms_av_oms2)
anova(coms_av_oms1,coms_av_oms3)

# prop germinated
coms_germ_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_germ_gms0,coms_germ_gms1)
anova(coms_germ_gms1,coms_germ_gms2)
anova(coms_germ_gms1,coms_germ_gms3)

coms_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_germ_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_germ_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
anova(coms_germ_oms0,coms_germ_oms1)
anova(coms_germ_oms1,coms_germ_oms2)
anova(coms_germ_oms1,coms_germ_oms3)

# seed-set
coms_ss_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_ss_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_ss_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
coms_ss_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS+import_nb_part_ID_out_co10   + (1|session:ID_full),family = binomial)
anova(coms_ss_gms0,coms_ss_gms1)
anova(coms_ss_gms1,coms_ss_gms2)
anova(coms_ss_gms1,coms_ss_gms3)

coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10* poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10+ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
coms_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor+import_nb_part_ID_out_co10    + (1|session:ID_full),family = binomial)
coms_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10    + (1|session:ID_full),family = binomial)
anova(coms_ss_oms0,coms_ss_oms1)
anova(coms_ss_oms1,coms_ss_oms2)
anova(coms_ss_oms1,coms_ss_oms3)

# seed weight
coms_weight_gms0 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms1 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms2 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms3 <- lme4::glmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+import_nb_part_ID_out_co10  + (1|session:ID_full))
anova(coms_weight_gms0,coms_weight_gms1)
anova(coms_weight_gms1,coms_weight_gms2)
anova(coms_weight_gms1,coms_weight_gms3)

coms_weight_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10   + (1|session:ID_full))
coms_weight_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(nb_part_ID_out_co10)), poids_sd ~ nb_part_ID_out_co10+import_nb_part_ID_out_co10   + (1|session:ID_full))
anova(coms_weight_oms0,coms_weight_oms1)
anova(coms_weight_oms1,coms_weight_oms2)
anova(coms_weight_oms1,coms_weight_oms3)
# 5 effects




# 5) flower scale as before but with individual data too --------

# table_for_gms <- dt %>% filter(gMS_fem_out >= 1) # only outcrossing
# table_for_gms <- dt %>% filter(nGenot > 1) # more than one seed genotyped
# table_for_gms <- dt %>% filter(!is.na(gMS_fem_out)) # all values
table_for_gms <- dt %>% filter((nGenot > 1)&gMS_fem_out >= 1) # more than one seed genotyped and outcrossing

#### prop avorted ####
# gms - nGenot as covariable
coms_av_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full),family = binomial)
coms_av_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + gMS_fem_out + nGenot  +  (1|session:ID_full), family = binomial)
coms_av_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + poll_treat_factor + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ poll_treat_factor + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
anova(coms_av_gms0,coms_av_gms1) # inter indiv * ttt
anova(coms_av_gms0,coms_av_gms2) # inter flower * ttt
anova(coms_av_gms3,coms_av_gms4) # ttt
anova(coms_av_gms3,coms_av_gms5) # indiv
anova(coms_av_gms3,coms_av_gms6) # flower

summary(coms_av_gms0)
car::Anova(coms_av_gms0) # 13 paramètres (11 fixe + 2 random) pour 222 observations
# flower : positive effect
performance::check_singularity(coms_av_gms3)
model_simple <- update(coms_av_gms3, . ~ . - (1|session))
anova(model_simple, coms_av_gms3, test = "Chisq")
# indeed removing session as random effect remove the warning
# intersession variance is null -> variance is driven by differences between individuals only
# we decide to keep it anyway (to be faithful to the data structure)

# oms - pl as covariable
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_av_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_av_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_av_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_av_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_av_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_av_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_av_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_av_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_av_oms0,coms_av_oms1) # inter indiv * ttt
anova(coms_av_oms0,coms_av_oms2) # inter flower * ttt
anova(coms_av_oms0,coms_av_oms3) # inter pl * ttt
anova(coms_av_oms4,coms_av_oms5) # ttt
anova(coms_av_oms4,coms_av_oms6) # indiv
anova(coms_av_oms4,coms_av_oms7) # flower
anova(coms_av_oms4,coms_av_oms8) # pl

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "high"))
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_av_oms0)
car::Anova(coms_av_oms0) # 13 paramètres (11 fixe + 2 random) pour 222 observations
# pl:
# low : negative
# medium : neutral
# high : neutral

performance::check_singularity(coms_av_oms0)
model_simple <- update(coms_av_oms0, . ~ . - (1|session))
anova(model_simple, coms_av_oms0, test = "Chisq")

#### seed-set ####
# gms - nGenot as covariable
coms_ss_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS * poll_treat_factor + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_ss_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full),family = binomial)
coms_ss_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS * poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_ss_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_ss_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_ss_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_ss_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
anova(coms_ss_gms0,coms_ss_gms1) # inter indiv * ttt
anova(coms_ss_gms0,coms_ss_gms2) # inter flower * ttt
anova(coms_ss_gms3,coms_ss_gms4) # ttt
anova(coms_ss_gms3,coms_ss_gms5) # indiv
anova(coms_ss_gms3,coms_ss_gms6) # flower

summary(coms_ss_gms3)
# ind : positive effect
# flow : negative effect

coms_av_oms0

# oms - pl as covariable
coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_ss_oms0,coms_ss_oms1) # inter indiv * ttt
anova(coms_ss_oms0,coms_ss_oms2) # inter flower * ttt
anova(coms_ss_oms0,coms_ss_oms3) # inter pl * ttt
anova(coms_ss_oms4,coms_ss_oms5) # ttt
anova(coms_ss_oms4,coms_ss_oms6) # indiv
anova(coms_ss_oms4,coms_ss_oms7) # flower
anova(coms_ss_oms4,coms_ss_oms8) # flower

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "high"))
coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_ss_oms0)
# ms flo:
# low : positive
# medium : positive
# high : neutral
# pl:
# low : negative
# medium : positive
# high : positive

#### prop germinated ####
# gms - nGenot as covariable
coms_germ_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS * poll_treat_factor + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_germ_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full),family = binomial)
coms_germ_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS * poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_germ_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_germ_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_germ_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_germ_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
anova(coms_germ_gms0,coms_germ_gms1) # inter indiv * ttt
anova(coms_germ_gms0,coms_germ_gms2) # inter flower * ttt
anova(coms_germ_gms3,coms_germ_gms4) # ttt
anova(coms_germ_gms3,coms_germ_gms5) # indiv
anova(coms_germ_gms3,coms_germ_gms6) # flower

# oms - pl as covariable
coms_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_germ_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_germ_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_germ_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_germ_oms0,coms_germ_oms1) # inter indiv * ttt
anova(coms_germ_oms0,coms_germ_oms2) # inter flower * ttt
anova(coms_germ_oms0,coms_germ_oms3) # inter pl * ttt
anova(coms_germ_oms4,coms_germ_oms5) # ttt
anova(coms_germ_oms4,coms_germ_oms6) # indiv
anova(coms_germ_oms4,coms_germ_oms7) # flower
anova(coms_germ_oms4,coms_germ_oms8) # pl

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "high"))
coms_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_germ_oms0)
# low : neutral (neg)
# medium : neutral (neg)
# high : neutral (pos)

#### seed weight ####
coms_weight_gms0 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms1 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms2 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms3 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+import_nb_part_ID_out_co10  + (1|session:ID_full))
anova(coms_weight_gms0,coms_weight_gms1)
anova(coms_weight_gms1,coms_weight_gms2)
anova(coms_weight_gms1,coms_weight_gms3)

coms_weight_oms0 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms1 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms2 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms3 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full))
coms_weight_oms4 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full))
coms_weight_oms5 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full))
coms_weight_oms6 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full))
coms_weight_oms7 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full))
coms_weight_oms8 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full))
anova(coms_weight_oms0,coms_weight_oms1) # inter indiv * ttt
anova(coms_weight_oms0,coms_weight_oms2) # inter flower * ttt
anova(coms_weight_oms0,coms_weight_oms3) # inter pl * ttt
anova(coms_weight_oms4,coms_weight_oms5) # ttt
anova(coms_weight_oms4,coms_weight_oms6) # indiv
anova(coms_weight_oms4,coms_weight_oms7) # flower
anova(coms_weight_oms4,coms_weight_oms8) # flower


# 6) same as bis but with ttt only as random effect and testing for it --------

# table_for_gms <- dt %>% filter(gMS_fem_out >= 1) # only outcrossing
# table_for_gms <- dt %>% filter(nGenot > 1) # more than one seed genotyped
# table_for_gms <- dt %>% filter(!is.na(gMS_fem_out)) # all values
table_for_gms <- dt %>% filter((nGenot > 1)&gMS_fem_out >= 1) # more than one seed genotyped and outcrossing

#### prop avorted ####
# gms - nGenot as covariable
coms_av_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms0rand <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out + nGenot + (1|poll_treat_factor) + (1|poll_treat_factor:session) + (1|poll_treat_factor:session:ID_full), family = binomial)
AIC(coms_av_gms0,coms_av_gms0rand)
# ok better as random effect

coms_av_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_fem_out + nGenot + (1|poll_treat_factor) + (1|poll_treat_factor:session) + (1|poll_treat_factor:session:ID_full),family = binomial)
coms_av_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + nGenot + (1|poll_treat_factor) + (1|poll_treat_factor:session) + (1|poll_treat_factor:session:ID_full), family = binomial)
coms_av_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out + (1|poll_treat_factor) + (1|poll_treat_factor:session) + (1|poll_treat_factor:session:ID_full), family = binomial)
anova(coms_av_gms0rand,coms_av_gms1) # indiv
anova(coms_av_gms0rand,coms_av_gms2) # fleur
anova(coms_av_gms0rand,coms_av_gms3) # ngenot

summary(coms_av_gms0rand)
# flower : positive effect

# oms - pl as covariable
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_av_oms0rand <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl + (1|poll_treat_factor) + (1|poll_treat_factor:session) + (1|poll_treat_factor:session:ID_full),family = binomial)
AIC(coms_av_oms0,coms_av_oms0rand)
# better with fixed effect

coms_av_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_av_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_av_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_av_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_av_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_av_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_av_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_av_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_av_oms0,coms_av_oms1) # inter indiv * ttt
anova(coms_av_oms0,coms_av_oms2) # inter flower * ttt
anova(coms_av_oms0,coms_av_oms3) # inter pl * ttt
anova(coms_av_oms4,coms_av_oms5) # ttt
anova(coms_av_oms4,coms_av_oms6) # indiv
anova(coms_av_oms4,coms_av_oms7) # flower
anova(coms_av_oms4,coms_av_oms8) # pl

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "high"))
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_av_oms0)
# pl:
# low : negative
# medium : neutral
# high : neutral

#### seed-set ####
# gms - nGenot as covariable
coms_ss_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS * poll_treat_factor + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_ss_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full),family = binomial)
coms_ss_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS * poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_ss_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_ss_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_ss_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_ss_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
anova(coms_ss_gms0,coms_ss_gms1) # inter indiv * ttt
anova(coms_ss_gms0,coms_ss_gms2) # inter flower * ttt
anova(coms_ss_gms3,coms_ss_gms4) # ttt
anova(coms_ss_gms3,coms_ss_gms5) # indiv
anova(coms_ss_gms3,coms_ss_gms6) # flower

summary(coms_ss_gms3)
# ind : positive effect
# flow : negative effect

# oms - pl as covariable
coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_ss_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_ss_oms0,coms_ss_oms1) # inter indiv * ttt
anova(coms_ss_oms0,coms_ss_oms2) # inter flower * ttt
anova(coms_ss_oms0,coms_ss_oms3) # inter pl * ttt
anova(coms_ss_oms4,coms_ss_oms5) # ttt
anova(coms_ss_oms4,coms_ss_oms6) # indiv
anova(coms_ss_oms4,coms_ss_oms7) # flower
anova(coms_ss_oms4,coms_ss_oms8) # flower

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "high"))
coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_ss_oms0)
# ms flo:
# low : positive
# medium : positive
# high : neutral
# pl:
# low : negative
# medium : positive
# high : positive

#### prop germinated ####
# gms - nGenot as covariable
coms_germ_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS * poll_treat_factor + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_germ_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + gMS_fem_out * poll_treat_factor + nGenot  + (1|session:ID_full),family = binomial)
coms_germ_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS * poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_germ_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + poll_treat_factor + gMS_fem_out + nGenot  + nGenot + (1|session:ID_full), family = binomial)
coms_germ_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_germ_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor + gMS_fem_out + nGenot  + (1|session:ID_full), family = binomial)
coms_germ_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
anova(coms_germ_gms0,coms_germ_gms1) # inter indiv * ttt
anova(coms_germ_gms0,coms_germ_gms2) # inter flower * ttt
anova(coms_germ_gms3,coms_germ_gms4) # ttt
anova(coms_germ_gms3,coms_germ_gms5) # indiv
anova(coms_germ_gms3,coms_germ_gms6) # flower

# oms - pl as covariable
coms_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_germ_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_germ_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
coms_germ_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full),family = binomial)
coms_germ_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full),family = binomial)
anova(coms_germ_oms0,coms_germ_oms1) # inter indiv * ttt
anova(coms_germ_oms0,coms_germ_oms2) # inter flower * ttt
anova(coms_germ_oms0,coms_germ_oms3) # inter pl * ttt
anova(coms_germ_oms4,coms_germ_oms5) # ttt
anova(coms_germ_oms4,coms_germ_oms6) # indiv
anova(coms_germ_oms4,coms_germ_oms7) # flower
anova(coms_germ_oms4,coms_germ_oms8) # pl

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "high"))
coms_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_germ_oms0)
# low : neutral (neg)
# medium : neutral (neg)
# high : neutral (pos)

#### seed weight ####
coms_weight_gms0 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS*poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms1 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms2 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ poll_treat_factor+import_nb_part_ID_out_co10  + (1|session:ID_full))
coms_weight_gms3 <- lme4::lmer(data = dt %>% filter(!is.na(gMS)), poids_sd ~ gMS+import_nb_part_ID_out_co10  + (1|session:ID_full))
anova(coms_weight_gms0,coms_weight_gms1)
anova(coms_weight_gms1,coms_weight_gms2)
anova(coms_weight_gms1,coms_weight_gms3)

coms_weight_oms0 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms1 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms2 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 + pl * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms3 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl  + (1|session:ID_full))
coms_weight_oms4 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full))
coms_weight_oms5 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full))
coms_weight_oms6 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ poll_treat_factor + import_nb_part_ID_out_co10 + pl  + (1|session:ID_full))
coms_weight_oms7 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + poll_treat_factor + pl  + (1|session:ID_full))
coms_weight_oms8 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), poids_sd ~ nb_part_ID_out_co10 + poll_treat_factor + import_nb_part_ID_out_co10  + (1|session:ID_full))
anova(coms_weight_oms0,coms_weight_oms1) # inter indiv * ttt
anova(coms_weight_oms0,coms_weight_oms2) # inter flower * ttt
anova(coms_weight_oms0,coms_weight_oms3) # inter pl * ttt
anova(coms_weight_oms4,coms_weight_oms5) # ttt
anova(coms_weight_oms4,coms_weight_oms6) # indiv
anova(coms_weight_oms4,coms_weight_oms7) # flower
anova(coms_weight_oms4,coms_weight_oms8) # flower



# 7) as analysis 5 but with pl in interaction for gms too --------

# table_for_gms <- dt %>% filter(gMS_fem_out >= 1) # only outcrossing
# table_for_gms <- dt %>% filter(nGenot > 1) # more than one seed genotyped
# table_for_gms <- dt %>% filter(!is.na(gMS_fem_out)) # all values
table_for_gms <- dt %>% 
  filter((nGenot > 1)&gMS_fem_out >= 1) %>% # more than one seed genotyped and outcrossing
  filter(!is.na(pl)) %>% # to keep the same dataset
  dplyr::mutate( # if convergence problem, center-scale continuous
    gMS_c = scale(gMS, center = TRUE, scale = TRUE),
    gMS_fem_out_c = scale(gMS_fem_out, center = TRUE, scale = TRUE),
    pl_c = scale(pl, center = TRUE, scale = TRUE),
    nGenot_c = scale(nGenot, center = TRUE, scale = TRUE)
  )
  

#### prop avorted ####
##### gms #####
# gms - nGenot as covariable

coms_av_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + gMS_fem_out * poll_treat_factor + pl * poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out * poll_treat_factor + pl * poll_treat_factor + nGenot  + (1|session:ID_full),family = binomial)
coms_av_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS * poll_treat_factor + gMS_fem_out + pl * poll_treat_factor + nGenot  +  (1|session:ID_full), family = binomial)
coms_av_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + poll_treat_factor + gMS_fem_out * poll_treat_factor + pl + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out + pl + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out + pl + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_fem_out + pl + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms7 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + pl + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms8 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out + poll_treat_factor + nGenot  + (1|session:ID_full), family = binomial)
coms_av_gms9 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS + gMS_fem_out + pl + poll_treat_factor  + (1|session:ID_full), family = binomial)

anova(coms_av_gms0,coms_av_gms1) # inter indiv * ttt
anova(coms_av_gms0,coms_av_gms2) # inter flower * ttt
anova(coms_av_gms0,coms_av_gms3) # inter pl * ttt
anova(coms_av_gms4,coms_av_gms5) # ttt
anova(coms_av_gms4,coms_av_gms6) # indiv
anova(coms_av_gms4,coms_av_gms7) # flower
anova(coms_av_gms4,coms_av_gms8) # pl
anova(coms_av_gms4,coms_av_gms9) # nGenot

summary(coms_av_gms0)
car::Anova(coms_av_gms0) # 14 paramètres (12 fixe + 2 random) pour 80 observations...
# flower : positive effect
performance::check_singularity(coms_av_gms0)
performance::check_collinearity(coms_av_gms0)
# okay il faut centrer-réduire vraiment, on regarde si ça enlève les warnings
# next try: 
coms_av_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c * poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full),family = binomial)
coms_av_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c * poll_treat_factor + gMS_fem_out_c + pl_c * poll_treat_factor + nGenot_c  +  (1|session:ID_full), family = binomial)
coms_av_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms7 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms8 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms9 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor  + (1|session:ID_full), family = binomial)

# again a problem, next try by removing session as random (variance =0):
# next try: 
coms_av_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c * poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full),family = binomial)
coms_av_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c * poll_treat_factor + gMS_fem_out_c + pl_c * poll_treat_factor + nGenot_c  +  (1|session:ID_full), family = binomial)
coms_av_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms7 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms8 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial)
coms_av_gms9 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor  + (1|session:ID_full), family = binomial)

# yes better, but removing session effect... mouaif
# try with another optimizer:
coms_av_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c * poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c * poll_treat_factor + gMS_fem_out_c + pl_c * poll_treat_factor + nGenot_c  +  (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms7 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms8 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_gms9 <- lme4::glmer(data = table_for_gms, cbind(nbGr, nbAv) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

performance::check_singularity(coms_av_gms0)
performance::check_collinearity(coms_av_gms0)
# no problem anymore

anova(coms_av_gms0,coms_av_gms1) # inter indiv * ttt
anova(coms_av_gms0,coms_av_gms2) # inter flower * ttt
anova(coms_av_gms0,coms_av_gms3) # inter pl * ttt
anova(coms_av_gms4,coms_av_gms5) # ttt
anova(coms_av_gms4,coms_av_gms6) # indiv
anova(coms_av_gms4,coms_av_gms7) # flower
anova(coms_av_gms4,coms_av_gms8) # pl
anova(coms_av_gms4,coms_av_gms9) # nGenot
# and not any effect

##### oms #####
# oms - pl as covariable 
# try another optimizer do not change anything
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c + poll_treat_factor + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_av_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
anova(coms_av_oms0,coms_av_oms1) # inter indiv * ttt
anova(coms_av_oms0,coms_av_oms2) # inter flower * ttt
anova(coms_av_oms0,coms_av_oms3) # inter pl * ttt
anova(coms_av_oms4,coms_av_oms5) # ttt
anova(coms_av_oms4,coms_av_oms6) # indiv
anova(coms_av_oms4,coms_av_oms7) # flower
anova(coms_av_oms4,coms_av_oms8) # pl

performance::check_singularity(coms_av_gms0)
performance::check_collinearity(coms_av_gms0)

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "medium"))
coms_av_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr,nbAv) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
summary(coms_av_oms0)
car::Anova(coms_av_oms0) # 13 paramètres (11 fixe + 2 random) pour 222 observations
# pl:
# low : negative
# medium : neutral
# high : neutral

# bon bilan : corriger les problèmes de convergence ne change pas drastiquement les pvalues ni les estimates
# on continue  en vérifiant qu'on a pas de problèmes de singularities,
# et qu'on ne descende pas en dessous de 80 points pour 14 paramètres comme vu précédemment

#### seed-set ####
##### gms #####
# gms - nGenot as covariable
coms_ss_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c * poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c * poll_treat_factor + gMS_fem_out_c + pl_c * poll_treat_factor + nGenot_c  +  (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c + poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c + gMS_fem_out_c + pl_c + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms7 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms8 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c + gMS_fem_out_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_gms9 <- lme4::glmer(data = table_for_gms, cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
performance::check_singularity(coms_ss_gms0) # ok
performance::check_collinearity(coms_ss_gms0) # ok

anova(coms_ss_gms0,coms_ss_gms1) # inter indiv * ttt
anova(coms_ss_gms0,coms_ss_gms2) # inter flower * ttt
anova(coms_ss_gms0,coms_ss_gms3) # inter pl * ttt
anova(coms_ss_gms4,coms_ss_gms5) # ttt
anova(coms_ss_gms4,coms_ss_gms6) # indiv
anova(coms_ss_gms4,coms_ss_gms7) # flower
anova(coms_ss_gms4,coms_ss_gms8) # pl
anova(coms_ss_gms4,coms_ss_gms9) # nGenot

summary(coms_ss_gms4) # 77 points

##### oms #####
# oms - pl as covariable
coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c + poll_treat_factor + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
anova(coms_ss_oms0,coms_ss_oms1) # inter indiv * ttt
anova(coms_ss_oms0,coms_ss_oms2) # inter flower * ttt
anova(coms_ss_oms0,coms_ss_oms3) # inter pl * ttt
anova(coms_ss_oms4,coms_ss_oms5) # ttt
anova(coms_ss_oms4,coms_ss_oms6) # indiv
anova(coms_ss_oms4,coms_ss_oms7) # flower
anova(coms_ss_oms4,coms_ss_oms8) # pl

performance::check_singularity(coms_ss_oms0) # ok
performance::check_collinearity(coms_ss_oms0) # ok

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "high"))
coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_SR,(nbOv - nbGr_SR)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_ss_oms0) # 182 points
# ms flo:
# low : positive
# medium : positive
# high : neutral
# pl:
# low : negative
# medium : positive
# high : positive

#### prop germinated ####
##### gms ##### 
# gms - nGenot as covariable
coms_germ_gms0 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c * poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms1 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms2 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c * poll_treat_factor + gMS_fem_out_c + pl_c * poll_treat_factor + nGenot_c  +  (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms3 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c + poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms4 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms5 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c + gMS_fem_out_c + pl_c + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms6 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms7 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms8 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c + gMS_fem_out_c + poll_treat_factor + nGenot_c  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_germ_gms9 <- lme4::glmer(data = table_for_gms, cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor  + (1|session:ID_full), family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
performance::check_singularity(coms_germ_gms0) # ok
performance::check_collinearity(coms_germ_gms0) # still high correlation for some variables

anova(coms_germ_gms0,coms_germ_gms1) # inter indiv * ttt
anova(coms_germ_gms0,coms_germ_gms2) # inter flower * ttt
anova(coms_germ_gms0,coms_germ_gms3) # inter pl * ttt
anova(coms_germ_gms4,coms_germ_gms5) # ttt
anova(coms_germ_gms4,coms_germ_gms6) # indiv
anova(coms_germ_gms4,coms_germ_gms7) # flower
anova(coms_germ_gms4,coms_germ_gms8) # pl
anova(coms_germ_gms4,coms_germ_gms9) # nGenot

##### oms #####
# oms - pl as covariable
coms_ss_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms1 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms2 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call * poll_treat_factor  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms3 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms4 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms5 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms6 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms7 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c + poll_treat_factor + pl_call  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
coms_ss_oms8 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c  + (1|session:ID_full),family = binomial,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
anova(coms_ss_oms0,coms_ss_oms1) # inter indiv * ttt
anova(coms_ss_oms0,coms_ss_oms2) # inter flower * ttt
anova(coms_ss_oms0,coms_ss_oms3) # inter pl * ttt
anova(coms_ss_oms4,coms_ss_oms5) # ttt
anova(coms_ss_oms4,coms_ss_oms6) # indiv
anova(coms_ss_oms4,coms_ss_oms7) # flower
anova(coms_ss_oms4,coms_ss_oms8) # pl

dt <- dt %>%
  mutate(poll_treat_factor=forcats::fct_relevel(poll_treat_factor, "low"))
coms_germ_oms0 <- lme4::glmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10)&!is.na(pl)), cbind(nbGr_germ,(nbGr_semis - nbGr_germ)) ~ nb_part_ID_out_co10 * poll_treat_factor + import_nb_part_ID_out_co10 * poll_treat_factor + pl * poll_treat_factor  + (1|session:ID_full),family = binomial)
summary(coms_germ_oms0)
# low : neutral (neg)
# medium : neutral (neg)
# high : neutral (pos)

#### seed weight ####
##### gms #####
coms_weight_gms0 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c * poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full))
coms_weight_gms1 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c + gMS_fem_out_c * poll_treat_factor + pl_c * poll_treat_factor + nGenot_c  + (1|session:ID_full))
coms_weight_gms2 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c * poll_treat_factor + gMS_fem_out_c + pl_c * poll_treat_factor + nGenot_c  +  (1|session:ID_full))
coms_weight_gms3 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c + poll_treat_factor + gMS_fem_out_c * poll_treat_factor + pl_c + nGenot_c  + (1|session:ID_full))
coms_weight_gms4 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full))
coms_weight_gms5 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c + gMS_fem_out_c + pl_c + nGenot_c  + (1|session:ID_full))
coms_weight_gms6 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_fem_out_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full))
coms_weight_gms7 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c + pl_c + poll_treat_factor + nGenot_c  + (1|session:ID_full))
coms_weight_gms8 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c + gMS_fem_out_c + poll_treat_factor + nGenot_c  + (1|session:ID_full))
coms_weight_gms9 <- lme4::lmer(data = table_for_gms, poids_sd ~ gMS_c + gMS_fem_out_c + pl_c + poll_treat_factor  + (1|session:ID_full))
performance::check_singularity(coms_weight_gms0) # ok
performance::check_collinearity(coms_weight_gms0) # ok

anova(coms_weight_gms0,coms_weight_gms1) # inter indiv * ttt
anova(coms_weight_gms0,coms_weight_gms2) # inter flower * ttt
anova(coms_weight_gms0,coms_weight_gms3) # inter pl * ttt
anova(coms_weight_gms4,coms_weight_gms5) # ttt
anova(coms_weight_gms4,coms_weight_gms6) # indiv
anova(coms_weight_gms4,coms_weight_gms7) # flower
anova(coms_weight_gms4,coms_weight_gms8) # pl
anova(coms_weight_gms4,coms_weight_gms9) # nGenot

summary(coms_ss_gms0) # 80 points

##### oms #####
# oms - pl as covariable
coms_weight_oms0 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms1 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms2 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call * poll_treat_factor  + (1|session:ID_full))
coms_weight_oms3 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c * poll_treat_factor + import_nb_part_ID_out_co10_c * poll_treat_factor + pl_call  + (1|session:ID_full))
coms_weight_oms4 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full))
coms_weight_oms5 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full))
coms_weight_oms6 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ poll_treat_factor + import_nb_part_ID_out_co10_c + pl_call  + (1|session:ID_full))
coms_weight_oms7 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c + poll_treat_factor + pl_call  + (1|session:ID_full))
coms_weight_oms8 <- lme4::lmer(data = dt %>% filter(!is.na(import_nb_part_ID_out_co10_c)&!is.na(pl_call)), poids_sd ~ nb_part_ID_out_co10_c + poll_treat_factor + import_nb_part_ID_out_co10_c  + (1|session:ID_full))
anova(coms_weight_oms0,coms_weight_oms1) # inter indiv * ttt
anova(coms_weight_oms0,coms_weight_oms2) # inter flower * ttt
anova(coms_weight_oms0,coms_weight_oms3) # inter pl * ttt
anova(coms_weight_oms4,coms_weight_oms5) # ttt
anova(coms_weight_oms4,coms_weight_oms6) # indiv
anova(coms_weight_oms4,coms_weight_oms7) # flower
anova(coms_weight_oms4,coms_weight_oms8) # pl

performance::check_singularity(coms_weight_oms0) # ok
performance::check_collinearity(coms_weight_oms0) # ok


# 8) final analyses -------------------------------------------------------


# plot id -----------------------------------------------------------------

dt_id_raw <- dt |>
  group_by(ID_full,poll_treat_factor) |>
  summarise(ms_gms = unique(gMS),
            ms_oms = unique(nb_part_ID_out_co10))

dt_id <- dt_id_raw |>
  pivot_longer(
    cols = starts_with("ms"),
    names_to = c("ms"),   
    names_pattern = "ms_(gms|oms)" 
  ) 

mu_id <- dt_id_raw |> 
  group_by(poll_treat_factor) |>
  summarise(mean_gms = mean(ms_gms,na.rm = T),
            sd_gms = sd(ms_gms,na.rm = T),
            mean_oms = mean(ms_oms,na.rm = T),
            sd_oms = sd(ms_oms,na.rm = T)) |>
  pivot_longer(
    cols = -poll_treat_factor,
    names_to = c(".value", "ms"),   
    names_pattern = "(mean|sd)_(.*)" 
  )

dt_id$ms <- factor(dt_id$ms, levels = c("oms", "gms"))

plot_id <- ggplot(data = dt_id, aes(x = value, fill = poll_treat_factor)) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(ms ~ ., ncol = 1, as.table = FALSE) +
  geom_histogram(aes(y=after_stat(density)),alpha=0.5) +
  geom_density(alpha=0.6) +
  ggtitle("Individual scale") +
  scale_fill_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  scale_color_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  geom_vline(data=mu_id, aes(xintercept=mean, color=poll_treat_factor),
           linetype="dashed")
plot_id

mu_id_oms <- dt_id_raw |> 
  group_by(poll_treat_factor) |>
  summarise(mean = mean(ms_oms,na.rm = T),
            sd = sd(ms_oms,na.rm = T)) 

plot_id_oms <- ggplot(data = dt_id %>% filter(ms == "oms"), aes(x = value, fill = poll_treat_factor)) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_histogram(aes(y=after_stat(density)),alpha=0.5) +
  geom_density(alpha=0.6) +
  ggtitle("Individual scale") +
  scale_fill_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  scale_color_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  geom_vline(data=mu_id_oms, aes(xintercept=mean, color=poll_treat_factor),
             linetype="dashed")
plot_id_oms


# plot flower -----------------------------------------------------------------

dt_flo_raw <- dt |>
  select(ID_full, poll_treat_factor, id_flow, gMS_fem_out, import_nb_part_ID_out_co10) |>
  rename(ms_gms = gMS_fem_out,
         ms_oms = import_nb_part_ID_out_co10)

dt_flo <- dt_flo_raw |>
  pivot_longer(
    cols = starts_with("ms"),
    names_to = c("ms"),   
    names_pattern = "ms_(gms|oms)" 
  ) 

mu_flo <- dt_flo_raw |> 
  group_by(poll_treat_factor) |>
  summarise(mean_gms = mean(ms_gms,na.rm = T),
            sd_gms = sd(ms_gms,na.rm = T),
            mean_oms = mean(ms_oms,na.rm = T),
            sd_oms = sd(ms_oms,na.rm = T)) |>
  pivot_longer(
    cols = -poll_treat_factor,
    names_to = c(".value", "ms"),   
    names_pattern = "(mean|sd)_(.*)" 
  )

dt_flo$ms <- factor(dt_flo$ms, levels = c("oms", "gms"))

plot_flo <- ggplot(data = dt_flo, aes(x = value, fill = poll_treat_factor)) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(ms ~ ., ncol = 1, as.table = FALSE) +
  geom_histogram(aes(y=after_stat(density)),alpha=0.5) +
  geom_density(alpha=0.6) +
  ggtitle("Flower scale") +
  scale_fill_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  scale_color_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  geom_vline(data=mu_flo, aes(xintercept=mean, color=poll_treat_factor),
             linetype="dashed")
plot_flo

mu_flo_oms <- dt_flo_raw |> 
  group_by(poll_treat_factor) |>
  summarise(mean = mean(ms_oms,na.rm = T),
            sd = sd(ms_oms,na.rm = T)) 

plot_flo_oms <- ggplot(data = dt_flo %>% filter(ms == "oms"), aes(x = value, fill = poll_treat_factor)) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_histogram(aes(y=after_stat(density)),alpha=0.5) +
  geom_density(alpha=0.6) +
  ggtitle("Flower scale") +
  scale_fill_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  scale_color_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  geom_vline(data=mu_flo_oms, aes(xintercept=mean, color=poll_treat_factor),
             linetype="dashed")
plot_flo_oms


# plot pollen load -----------------------------------------------------------------

dt_pl <- dt |>
  select(ID_full, poll_treat_factor, id_flow, pl) 

mu_pl <- dt_pl |> 
  group_by(poll_treat_factor) |>
  summarise(mean = mean(pl,na.rm = T),
            sd = sd(pl,na.rm = T)) 

plot_pl <- ggplot(data = dt_pl, aes(x = pl, fill = poll_treat_factor)) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_histogram(aes(y=after_stat(density)),alpha=0.5) +
  geom_density(alpha=0.6) +
  ggtitle("Flower scale - Pollen load") +
  scale_fill_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  scale_color_manual(values = c("#A3333D","#CA755E","#f0b67f")) +
  geom_vline(data=mu_pl, aes(xintercept=mean, color=poll_treat_factor),
             linetype="dashed")
plot_pl
