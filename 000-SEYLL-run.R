### BeBOD YLL // RUN ALL YEARS // 2019
### 03/08/2022

## settings
n <- 100L
options(future.globals.maxSize = 1024 ^ 3)

## load helpers
setwd("//sciensano.be/fs/1140_DATA/DALY/08_BeBOD/04_YLL/04_ANALYSES/2019")
source("00-SEYLL-helpers.R")
save.image("settings.RData")

## load full dataset
rmarkdown::render("01-SEYLL-create.R")

## save, clean, reload
saveRDS(mrt_all, file = "mrt2019.rds")
rm(list = ls())
load("settings.RData"); unlink("settings.RData")
mrt_all <- readRDS("mrt2019.rds")

## run simulations for all years
all_yrs <- 2004:2019
for (i in seq_along(all_yrs)) {
  yrs <- seq(all_yrs[i] - 4, all_yrs[i])
  rmarkdown::render(
    "02-SEYLL-redistribute.R",
    output_file = sprintf("02-SEYLL-redistribute-%s.html", tail(yrs, 1)))
}

## TO DO // define U codes as _gc

## compile results
rm(mrt_all)
source("03-SEYLL-results.R")
source("04-SEYLL-compile.R")

## make plots
#source("05-SEYLL-plots.R")