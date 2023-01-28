#' ---
#' title: BeBOD SEYLL `r tail(yrs, 1)`
#' subtitle: MCOD redistribution model
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#' ---

#' # Settings

## make sure we can reproduce these results
set.seed(264)

## subset complete dataset
yrs
mrt <- subset(mrt_all, YEAR %fin% yrs)


#' # Export IDDs

mrt_last <- subset(mrt, YEAR == tail(yrs, 1))
is_idd <- map_gbd(mrt_last$UCAUSE4_2)[, "cause1"] == "Garbage code"
idd_tab <- table(mrt_last$UCAUSE4_2[is_idd])
idd_tab <- as.data.frame(idd_tab)
idd_tab$ICD_LAB <- idd$ICD_LAB[match(idd_tab$Var1, idd$ICD_CODE)]
names(idd_tab) <- c("ICD_CODE", "N", "ICD_LAB")
idd_tab <- idd_tab[, c("ICD_CODE", "ICD_LAB", "N")]
idd_tab$ICD_CODE <- as.character(idd_tab$ICD_CODE)

saveWidget(
  datatable(
    idd_tab,
    rownames = FALSE,
    filter = "top",
    options = list(pageLength = nrow(idd_tab)),
    caption = htmltools::tags$caption(
      style = "text-align: left;",
      sprintf("Ill-defined deaths for the year %s", tail(yrs, 1)))),
  file = sprintf("idd-%s.html", tail(yrs, 1)))


#' # Prepare redistributions

#############################################################################
## STEP 0 > TAKE SPECIFIC CODES #############################################
#############################################################################

## prepare redistribution dataframe
red <-
  mrt[, c("CLE", "YEAR", "AGE", "SEX", "AGEGRPSEX", "REGIOJ", "PROV", "GSC",
          "YLL", "UCAUSE4_2")]

## assign specific codes
red$ICD_RED0 <- NA
red$ICD_RED0[mrt$CAUSE1 != "Garbage code"] <-
  mrt$UCAUSE4_2[mrt$CAUSE1 != "Garbage code"]

## assign drug-related deaths cf EMCDDA
red$ICD_RED0[mrt$MIIICVC2 == "X41" & mrt$MIIICVC4D == "T436"] <- "F19"
red$ICD_RED0[mrt$MIIICVC2 == "X42" & mrt$MIIICVC4D == "T40"] <- "F19"
red$ICD_RED0[mrt$MIIICVC2 == "X44" & mrt$MIIICVC4D == "T436"] <- "F19"
red$ICD_RED0[mrt$MIIICVC2 == "X61" & mrt$MIIICVC4D == "T436"] <- "F19"
red$ICD_RED0[mrt$MIIICVC2 == "X62" & mrt$MIIICVC4D == "T40"] <- "F19"
red$ICD_RED0[mrt$MIIICVC2 == "X64" & mrt$MIIICVC4D == "T436"] <- "F19"
red$ICD_RED0[mrt$MIIICVC2 == "Y11" & mrt$MIIICVC4D == "T436"] <- "F19"
red$ICD_RED0[mrt$MIIICVC2 == "Y12" & mrt$MIIICVC4D == "T40"] <- "F19"
sum(red$ICD_RED0 == "F19", na.rm = TRUE) - sum(red$UCAUSE4_2 == "F19")

## summary
sum(!is.na(red$ICD_RED0))
sum(!is.na(red$ICD_RED0)) / nrow(red)

## identify codes to exclude from redistribution process
red$GBD3_RED0 <- map_gbd(red$ICD_RED0)[, "cause3"]
gbd3_tab <- table(red$GBD3_RED0)
(gbd3_exclude <- names(gbd3_tab[gbd3_tab < 5]))  # less than 5 (1/yr)
red$INCLUDE <- !(red$GBD3_RED0 %fin% gbd3_exclude)
sum(!red$INCLUDE); sum(gbd3_tab[gbd3_exclude]) # same?
red$GBD3_RED0 <- NULL

#' # Define sim function

sim_mrt <-
function() {

#############################################################################
## STEP 1a > ICD GROUPS #####################################################
#############################################################################

red$ICD_RED1 <- NA
red$ICD_RED1[!is.na(red$ICD_RED0)] <- red$ICD_RED0[!is.na(red$ICD_RED0)]

icd <- na.omit(unique(idd2$TARGET[idd2$REDISTRIBUTION == "ICD"]))

#i <- 2
#i <- 56  # stroke

for (i in seq_along(icd)) {

  # extract ill-defined deaths
  iddi <- subset(idd2, TARGET == icd[i])
  redi_id <- red$UCAUSE4_2 %fin% iddi[, "ICD_CODE"] & is.na(red$ICD_RED1)
  head(red[redi_id, ])

  # extract target
  targetICD <- explode_icd(icd[i])
  target <-
    subset(red, INCLUDE & grepl(targetICD, red$ICD_RED1, perl = TRUE))
  target_list <- with(target, tapply(ICD_RED1, AGEGRPSEX, list))

  # if target empty, use all ages
  target_list_mis_m <-
    which(grepl("M", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_m] <-
    list(na.omit(unlist(unname(target_list[grepl("M", names(target_list))]))))
  target_list_mis_f <-
    which(grepl("F", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_f] <-
    list(na.omit(unlist(unname(target_list[grepl("F", names(target_list))]))))

  # replace by target
  for (agesex in levels(red$AGEGRPSEX)) {
    redi_id_as <- redi_id & red$AGEGRPSEX == agesex
    if (sum(redi_id_as) > 0)
      red[redi_id_as, "ICD_RED1"] <-
        dqsample(target_list[[agesex]], sum(redi_id_as), replace = TRUE)
  }
}

tail(sort(table(red$ICD_RED0)), 10)
tail(sort(table(red$ICD_RED1)), 10)

sum(!is.na(red$ICD_RED1)) - sum(!is.na(red$ICD_RED0))
(sum(!is.na(red$ICD_RED1)) - sum(!is.na(red$ICD_RED0))) / nrow(red)


#############################################################################
## STEP 1b > GBD GROUPS #####################################################
#############################################################################

### LEVEL 4 ----------------------------------------------------------------

L4 <- na.omit(unique(idd2$TARGET[idd2$REDISTRIBUTION == "Level4"]))

for (i in seq_along(L4)) {

  # update GBD codes
  red$GBD_RED1 <-
    causelist$Code[match(map_gbd(red$ICD_RED1)[, 1], causelist$Level4)]

  # extract ill-defined deaths
  iddi <- subset(idd2, TARGET == L4[i])
  redi_id <- red$UCAUSE4_2 %fin% iddi[, "ICD_CODE"] & is.na(red$ICD_RED1)
  head(red[redi_id, ])

  # extract target
  targetICD <- gsub(" ", "", unlist(strsplit(L4[i], "\\|")))
  targetICD <- paste(targetICD , collapse = "|")
  target <-
    subset(red, INCLUDE & grepl(targetICD, red$GBD_RED1, perl = TRUE))
  target_list <- with(target, tapply(ICD_RED1, AGEGRPSEX, list))

  # if target empty, use all ages
  target_list_mis_m <-
    which(grepl("M", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_m] <-
    list(na.omit(unlist(unname(target_list[grepl("M", names(target_list))]))))
  target_list_mis_f <-
    which(grepl("F", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_f] <-
    list(na.omit(unlist(unname(target_list[grepl("F", names(target_list))]))))

  # replace by target
  for (agesex in levels(red$AGEGRPSEX)) {
    redi_id_as <- redi_id & red$AGEGRPSEX == agesex
    if (sum(redi_id_as) > 0)
      red[redi_id_as, "ICD_RED1"] <-
        dqsample(target_list[[agesex]], sum(redi_id_as), replace = TRUE)
  }
}

### LEVEL 3 ----------------------------------------------------------------

L3 <- na.omit(unique(idd2$TARGET[idd2$REDISTRIBUTION == "Level3"]))

for (i in seq_along(L3)) {

  # update GBD codes
  red$GBD_RED1 <-
    causelist$Code[match(map_gbd(red$ICD_RED1)[, 1], causelist$Level4)]

  # extract ill-defined deaths
  iddi <- subset(idd2, TARGET == L3[i])
  redi_id <- red$UCAUSE4_2 %fin% iddi[, "ICD_CODE"] & is.na(red$ICD_RED1)
  head(red[redi_id, ])

  # extract target
  targetICD <- gsub(" ", "", unlist(strsplit(L3[i], "\\|")))
  targetICD <- paste(targetICD , collapse = "|")
  target <-
    subset(red, INCLUDE & grepl(targetICD, red$GBD_RED1, perl = TRUE))
  target_list <- with(target, tapply(ICD_RED1, AGEGRPSEX, list))

  # if target empty, use all ages
  target_list_mis_m <-
    which(grepl("M", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_m] <-
    list(na.omit(unlist(unname(target_list[grepl("M", names(target_list))]))))
  target_list_mis_f <-
    which(grepl("F", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_f] <-
    list(na.omit(unlist(unname(target_list[grepl("F", names(target_list))]))))

  # replace by target
  for (agesex in levels(red$AGEGRPSEX)) {
    redi_id_as <- redi_id & red$AGEGRPSEX == agesex
    if (sum(redi_id_as) > 0)
      red[redi_id_as, "ICD_RED1"] <-
        dqsample(target_list[[agesex]], sum(redi_id_as), replace = TRUE)
  }
}

### LEVEL 2 ----------------------------------------------------------------

L2 <- na.omit(unique(idd2$TARGET[idd2$REDISTRIBUTION == "Level2"]))

for (i in seq_along(L2)) {

  # update GBD codes
  red$GBD_RED1 <-
    causelist$Code[match(map_gbd(red$ICD_RED1)[, 1], causelist$Level4)]

  # extract ill-defined deaths
  iddi <- subset(idd2, TARGET == L2[i])
  redi_id <- red$UCAUSE4_2 %fin% iddi[, "ICD_CODE"] & is.na(red$ICD_RED1)
  head(red[redi_id, ])

  # extract target
  targetICD <- gsub(" ", "", unlist(strsplit(L2[i], "\\|")))
  targetICD <- paste(targetICD , collapse = "|")
  target <-
    subset(red, INCLUDE & grepl(targetICD, red$GBD_RED1, perl = TRUE))
  target_list <- with(target, tapply(ICD_RED1, AGEGRPSEX, list))

  # if target empty, use all ages
  target_list_mis_m <-
    which(grepl("M", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_m] <-
    list(na.omit(unlist(unname(target_list[grepl("M", names(target_list))]))))
  target_list_mis_f <-
    which(grepl("F", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_f] <-
    list(na.omit(unlist(unname(target_list[grepl("F", names(target_list))]))))

  # replace by target
  for (agesex in levels(red$AGEGRPSEX)) {
    redi_id_as <- redi_id & red$AGEGRPSEX == agesex
    if (sum(redi_id_as) > 0)
      red[redi_id_as, "ICD_RED1"] <-
        dqsample(target_list[[agesex]], sum(redi_id_as), replace = TRUE)
  }
}

red$GBD_RED1 <- NULL

tail(sort(table(red$ICD_RED0)), 10)
tail(sort(table(red$ICD_RED1)), 10)

sum(!is.na(red$ICD_RED1)) - sum(!is.na(red$ICD_RED0))
(sum(!is.na(red$ICD_RED1)) - sum(!is.na(red$ICD_RED0))) / nrow(red)


#############################################################################
## STEP 2 > PACKAGES ########################################################
#############################################################################

red$ICD_RED2 <- NA
red$ICD_RED2[!is.na(red$ICD_RED1)] <- red$ICD_RED1[!is.na(red$ICD_RED1)]

pkg <- na.omit(unique(idd$PKG))

for (i in seq_along(pkg)) {

  # extract ill-defined deaths
  iddi <- subset(idd, PKG == pkg[i])
  redi_id <- red$UCAUSE4_2 %fin% iddi[, "ICD_CODE"] & is.na(red$ICD_RED2)
  head(red[redi_id, ])

  # extract target > underlying where iddi is intermediate
  targetICD_id <-
    (mrt$MCOD_ICD_01 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_02 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_03 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_04 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_05 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_06 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_07 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_08 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_09 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_10 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_11 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_12 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_13 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_14 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_15 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_16 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_17 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_18 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_19 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_20 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_21 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_22 %fin% iddi$ICD_CODE |
     mrt$MCOD_ICD_23 %fin% iddi$ICD_CODE) &
    !(mrt$UCAUSE4 %fin% iddi$ICD_CODE) &
    !is.na(red$ICD_RED2); sum(targetICD_id)
  target <- subset(red, INCLUDE & targetICD_id)
  target_list <- with(target, tapply(ICD_RED2, AGEGRPSEX, list))

  # if target totally empty, use all causes
  if (nrow(target) == 0) {
    target <- subset(red, INCLUDE & !is.na(ICD_RED2))
    target_list <- with(target, tapply(ICD_RED2, AGEGRPSEX, list))
  }

  # if stratum-specific target empty, use all ages
  target_list_mis_m <-
    which(grepl("M", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_m] <-
    list(na.omit(unlist(unname(target_list[grepl("M", names(target_list))]))))
  target_list_mis_f <-
    which(grepl("F", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_f] <-
    list(na.omit(unlist(unname(target_list[grepl("F", names(target_list))]))))
  #str(target_list)

  # if still empty, use all sexes and ages
  target_list_mis_m <-
    which(grepl("M", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_m] <-
    list(na.omit(unlist(unname(target_list))))
  target_list_mis_f <-
    which(grepl("F", names(target_list)) & sapply(target_list, is.logical))
  target_list[target_list_mis_f] <-
    list(na.omit(unlist(unname(target_list))))
  #str(target_list)

  # remove sex-specific codes
  icd_m <- c("C61", "C610", "C619", "D291", "D400") # prostate cancer
  target_list_f <- which(grepl("F", names(target_list)))
  target_list[target_list_f] <-
    lapply(target_list[target_list_f], function(x) x[!(x %in% icd_m)])

  # replace by target
  for (agesex in levels(red$AGEGRPSEX)) {
    redi_id_as <- redi_id & red$AGEGRPSEX == agesex
    if (sum(redi_id_as) > 0)
      red[redi_id_as, "ICD_RED2"] <-
        dqsample(target_list[[agesex]], sum(redi_id_as), replace = TRUE)
  }
}

sum(!is.na(red$ICD_RED2)) - sum(!is.na(red$ICD_RED1))
(sum(!is.na(red$ICD_RED2)) - sum(!is.na(red$ICD_RED1))) / nrow(red)


#############################################################################
## STEP 3 > INTERNAL REDISTRIBUTION #########################################
#############################################################################

## identify deaths with specific code in intermediate ICD
id <-
  mrt$CAUSE4 == "Garbage code" &
  sapply(mrt$MCOD_GBD, function(x) !all(x == "_gc")) &
  is.na(red$ICD_RED2)
sort(table(unlist(mrt[id, "MCOD_GBD"])))

sum(id)                                     # nr of deaths
sum(id) / nrow(mrt)                         # prop of total deaths
sum(id) / sum(mrt$CAUSE4 == "Garbage code") # prop of ill-defined deaths

## per death, randomly assign one of the specific codes

red$ICD_RED3 <- NA
red$ICD_RED3[!is.na(red$ICD_RED2)] <- red$ICD_RED2[!is.na(red$ICD_RED2)]

red[id, ]$ICD_RED3 <- apply(mrt[id, ], 1, sim_internal)

tail(sort(table(red$ICD_RED2)), 10)
tail(sort(table(red$ICD_RED3)), 10)

sum(!is.na(red$ICD_RED3)) - sum(!is.na(red$ICD_RED2))
(sum(!is.na(red$ICD_RED3)) - sum(!is.na(red$ICD_RED2))) / nrow(red)

## update codes to exclude from redistribution process

red$GBD3_RED3 <- map_gbd(red$ICD_RED3)[, "cause3"]
gbd3_tab <- table(red$GBD3_RED3)
(gbd3_exclude <- names(gbd3_tab[gbd3_tab < 5]))  # less than 5 (1/yr)
red$INCLUDE <- !(red$GBD3_RED3 %fin% gbd3_exclude)
sum(!red$INCLUDE); sum(gbd3_tab[gbd3_exclude]) # same?

red$GBD3_RED3 <- NULL


#############################################################################
## STEP 4 > ALL CODES #######################################################
#############################################################################

red$ICD_RED4 <- NA
red$ICD_RED4[!is.na(red$ICD_RED3)] <- red$ICD_RED3[!is.na(red$ICD_RED3)]

L0 <- na.omit(unique(idd2$TARGET[idd$REDISTRIBUTION == "Level0"]))

# extract ill-defined deaths
iddi <- subset(idd2, TARGET == "ALL")
redi_id <- red$UCAUSE4_2 %fin% iddi[, "ICD_CODE"] & is.na(red$ICD_RED4)
head(red[redi_id, ])

# extract target
target <- subset(red, INCLUDE & !is.na(ICD_RED4))
target_list <- with(target, tapply(ICD_RED4, AGEGRPSEX, list))

# replace by target
for (agesex in levels(red$AGEGRPSEX)) {
  redi_id_as <- redi_id & red$AGEGRPSEX == agesex
  if (sum(redi_id_as) > 0)
    red[redi_id_as, "ICD_RED4"] <-
      dqsample(target_list[[agesex]], sum(redi_id_as), replace = TRUE)
}

tail(sort(table(red$ICD_RED0)), 12)
tail(sort(table(red$ICD_RED1)), 12)
tail(sort(table(red$ICD_RED2)), 12)
tail(sort(table(red$ICD_RED3)), 12)
tail(sort(table(red$ICD_RED4)), 12)

sum(!is.na(red$ICD_RED4)) - sum(!is.na(red$ICD_RED3))
(sum(!is.na(red$ICD_RED4)) - sum(!is.na(red$ICD_RED3))) / nrow(red)

sum(is.na(red$ICD_RED4))            # should be 0
head(subset(red, is.na(ICD_RED4)))  # should be empty

#############################################################################

## only keep last year
red_last <- subset(red, YEAR == max(red$YEAR))
nrow(red_last)

## map to ICD
red_last$GBD3_RED0 <-
  factor(map_gbd(red_last$ICD_RED0)[, "cause3"], unique(causelist$Level3))
red_last$GBD3_RED1 <-
  factor(map_gbd(red_last$ICD_RED1)[, "cause3"], unique(causelist$Level3))
red_last$GBD3_RED2 <-
  factor(map_gbd(red_last$ICD_RED2)[, "cause3"], unique(causelist$Level3))
red_last$GBD3_RED3 <-
  factor(map_gbd(red_last$ICD_RED3)[, "cause3"], unique(causelist$Level3))
red_last$GBD3_RED4 <-
  factor(map_gbd(red_last$ICD_RED4)[, "cause3"], unique(causelist$Level3))
red_last$GBD2_RED4 <-
  factor(map_gbd(red_last$ICD_RED4)[, "cause2"], unique(causelist$Level2))
red_last$GBD1_RED4 <-
  factor(map_gbd(red_last$ICD_RED4)[, "cause1"], unique(causelist$Level1))

## check if any NA codes
is_na <- is.na(red_last$GBD3_RED4)
if (any(is_na)) {
  saveRDS(subset(red_last, is_na),
          file = sprintf("01_SIM/mrt-sim-%s-na-%s.rds",
                         tail(yrs, 1),
                         format(Sys.time(), "%Y%m%d%H%M%S")))
  warning("Missing codes found in final dataset.")
}

## check if any garbage codes
is_idd <- red_last$GBD3_RED4 == "Garbage code"
if (any(is_idd)) {
  saveRDS(subset(red_last, is_idd),
          file = sprintf("01_SIM/mrt-sim-%s-gc-%s.rds",
                         tail(yrs, 1),
                         format(Sys.time(), "%Y%m%d%H%M%S")))
  warning("Garbage codes found in final dataset.")
}

return(red_last)
}


#' # Run simulations

plan(multisession)
mrt_sim <- future_replicate(n, sim_mrt(), simplify = FALSE)
saveRDS(mrt_sim, file = sprintf("01_SIM/mrt-sim-%s.rds", tail(yrs, 1)))
future:::ClusterRegistry("stop")


#' # Results

tab0 <- table(mrt_sim[[1]]$GBD3_RED0)
tab1 <- table(subset(mrt_sim[[1]], is.na(GBD3_RED0))$GBD3_RED1)
tab2 <- table(subset(mrt_sim[[1]], is.na(GBD3_RED1))$GBD3_RED2)
tab3 <- table(subset(mrt_sim[[1]], is.na(GBD3_RED2))$GBD3_RED3)
tab4 <- table(subset(mrt_sim[[1]], is.na(GBD3_RED3))$GBD3_RED4)

df <-
data.frame(
  CAUSE = c(names(tab0), names(tab1), names(tab2), names(tab3), names(tab4)),
  STEP = rep(0:4, times = c(length(tab0), length(tab1), length(tab2), length(tab3), length(tab4))),
  VALUE = c(tab0, tab1, tab2, tab3, tab4))
df$STEP <- factor(df$STEP, rev(unique(df$STEP)))

df_agg <- aggregate(VALUE ~ CAUSE, df, sum)
df_agg <- df_agg[order(df_agg$VALUE), ]

df$CAUSE <- factor(df$CAUSE, unique(df_agg$CAUSE))

df <- subset(df, CAUSE %fin% tail(df_agg$CAUSE, 50))


#' ## Redistributions
red_count <-
sapply(mrt_sim,
  function(x)
    c(sum(!is.na(x$ICD_RED0)),
      sum(!is.na(x$ICD_RED1)) - sum(!is.na(x$ICD_RED0)),
      sum(!is.na(x$ICD_RED2)) - sum(!is.na(x$ICD_RED1)),
      sum(!is.na(x$ICD_RED3)) - sum(!is.na(x$ICD_RED2)),
      sum(!is.na(x$ICD_RED4)) - sum(!is.na(x$ICD_RED3))))
red_count_tab <-
  cbind(TOTAL = rowMeans(red_count),
        PCT = rowMeans(t(t(red_count) / colSums(red_count))))
rownames(red_count_tab) <- c("SPEC", "ICD/GBD", "PKG", "INT", "ALL")
knitr::kable(red_count_tab)


#' ## Plot
#+ fig.height=12, fig.width=12
ggplot(df, aes(x = CAUSE, y = VALUE, fill = STEP)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_brewer(
    palette = "Spectral",
    labels = rev(c("Specific", "ICD/GBD", "PKG", "Internal", "ALL"))) +
  scale_y_continuous(NULL, labels = scales::comma) +
  scale_x_discrete(NULL)


#' ## Table
knitr::kable(row.names = FALSE, summarize(~ GBD3_RED4, sort = TRUE))
knitr::kable(row.names = FALSE, summarize(~ GBD2_RED4, sort = TRUE))


#' ## Table YLL
knitr::kable(row.names = FALSE, summarize(YLL ~ GBD3_RED4, sort = TRUE))
knitr::kable(row.names = FALSE, summarize(YLL ~ GBD2_RED4, sort = TRUE))


##rmarkdown::render("02-SEYLL-redistribute.R")