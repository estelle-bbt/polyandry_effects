targets::tar_load(results_plot)
results_plot$plot_q_oms
results_plot$plot_wo_contact
results_plot$plot_wo_contact_pl
results_plot$plot_w_contact
results_plot$plot_w_contact_pl
results_plot$plot_w_contact_contact

# wo contact
targets::tar_load(stat_flower_wo_contact)
View(stat_flower_wo_contact$estimate_table)

targets::tar_load(stat_id_wo_contact)
View(stat_id_wo_contact$estimate_table)

# w contact
targets::tar_load(stat_flower_w_contact)
View(stat_flower_w_contact$estimate_table)

targets::tar_load(stat_id_w_contact)
View(stat_id_w_contact$estimate_table)
