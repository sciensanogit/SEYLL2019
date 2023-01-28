### BeBOD YLL // HELPERS
### 20/08/2022

## required packages
library(arrow)
library(bd)
library(BeBOD)
library(DT)
library(dqrng)
library(fastmatch)
library(future.apply)
library(ggplot2)
library(htmltools)
library(htmlwidgets)
library(knitr)

## helper data

## import ill-defined code definitions
f <- "../../03_INPUTS/BeBOD-IDD-20220731.xlsx"
idd <- readxl(f)
idd2 <- subset(idd, is.na(PKG))

## import GBD cause list
gbd <- BeBOD:::gbd
causelist <- BeBOD:::causelist

## helper functions

## internal redistribution
sim_internal <-
function(x) {
  j <- x["MCOD_GBD"][[1]] != "_gc"
  sample(x["MCOD_ICD"][[1]][j], 1)
}

## summarize 'mrt'_sim list
summarize <-
function(f, sort = FALSE) {
  out_lst <-
  lapply(mrt_sim,
    function(x) as.data.frame(xtabs(f, x)))
  out_lab <- out_lst[[1]]
  out_lab$Freq <- NULL
  out_val <- sapply(out_lst, function(x) x$Freq)
  out_val <- t(apply(out_val, 1, mean_ci))
  colnames(out_val) <- c("VALUE", "VALUE_LWR", "VALUE_UPR")
  out <- cbind(out_lab, out_val)
  if (sort) out <- out[order(out$VALUE, decreasing = TRUE), ]
  return(out)
}