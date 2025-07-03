#' Project builder

## Install Dependencies (listed in DESCRIPTION) ----

remotes::install_deps(upgrade = "never")

## Run Project ----

targets::tar_make()


