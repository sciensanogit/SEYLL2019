### BEBOD / COMPILE SEYLL RESULTS

###
### LEVEL 0
###

## import results
mrt0 <-
lapply(all_yrs,
  function(y) readRDS(sprintf("02_RDS/mrt-yll-level0-%s.rds", y)))
mrt0 <- do.call("rbind", mrt0)

## save combined results
saveRDS(
  mrt0,
  file = sprintf("02_RDS/mrt-yll-level0-%s-%s.rds",
                 head(all_yrs, 1), tail(all_yrs, 1)))


###
### LEVEL 1
###

## import results
mrt1 <-
lapply(all_yrs,
  function(y) readRDS(sprintf("02_RDS/mrt-yll-level1-%s.rds", y)))
mrt1 <- do.call("rbind", mrt1)

## save combined results
saveRDS(
  mrt1,
  file = sprintf("02_RDS/mrt-yll-level1-%s-%s.rds",
                 head(all_yrs, 1), tail(all_yrs, 1)))


###
### LEVEL 2
###

## import results
mrt2 <-
lapply(all_yrs,
  function(y) readRDS(sprintf("02_RDS/mrt-yll-level2-%s.rds", y)))
mrt2 <- do.call("rbind", mrt2)

## save combined results
saveRDS(
  mrt2,
  file = sprintf("02_RDS/mrt-yll-level2-%s-%s.rds",
                 head(all_yrs, 1), tail(all_yrs, 1)))


###
### LEVEL 3
###

## import results
mrt3 <-
lapply(all_yrs,
  function(y) readRDS(sprintf("02_RDS/mrt-yll-level3-%s.rds", y)))
mrt3 <- do.call("rbind", mrt3)

## save combined results
saveRDS(
  mrt3,
  file = sprintf("02_RDS/mrt-yll-level3-%s-%s.rds",
                 head(all_yrs, 1), tail(all_yrs, 1)))

###
### ALL LEVELS
###

mrt0 <- cbind(LEVEL = 0, CAUSE = "ALL CAUSES", mrt0)
str(mrt0)

mrt1 <- cbind(LEVEL = 1, mrt1)
names(mrt1)[names(mrt1) == "CAUSE1"] <- "CAUSE"
str(mrt1)

mrt2$CAUSE1 <- NULL
mrt2 <- cbind(LEVEL = 2, mrt2)
names(mrt2)[names(mrt2) == "CAUSE2"] <- "CAUSE"
str(mrt2)

mrt3$CAUSE2 <- NULL
mrt3$CAUSE1 <- NULL
mrt3 <- cbind(LEVEL = 3, mrt3)
names(mrt3)[names(mrt3) == "CAUSE3"] <- "CAUSE"
str(mrt3)

mrt <- rbind(mrt0, mrt1, mrt2, mrt3)

## add PARENT cause
mrt$PARENT <- NA
mrt$PARENT[mrt$LEVEL == 1] <- "ALL CAUSES"
mrt$PARENT[mrt$LEVEL == 2] <-
  BeBOD:::causelist$Level1[
    match(mrt$CAUSE[mrt$LEVEL == 2], BeBOD:::causelist$Level2)]
mrt$PARENT[mrt$LEVEL == 3] <-
  BeBOD:::causelist$Level2[
    match(mrt$CAUSE[mrt$LEVEL == 3], BeBOD:::causelist$Level3)]
mrt <-
  mrt[, c("LEVEL", "CAUSE", "PARENT", "YEAR", "AGEGRP", "SEX", "REGIOJ",
          "MEASURE", "METRIC", "VAL_MEAN", "VAL_LWR", "VAL_UPR")]

## relevel factors
mrt$AGEGRP <-
  factor(mrt$AGEGRP,
         levels = c("[0,5)",   "[5,10)",  "[10,15)", "[15,20)", "[20,25)",
                    "[25,30)", "[30,35)", "[35,40)", "[40,45)", "[45,50)",
                    "[50,55)", "[55,60)", "[60,65)", "[65,70)", "[70,75)",
                    "[75,80)", "[80,85)", "[85,Inf)", "ALL", "BSP", "ESP"))

## save combined results
saveRDS(
  mrt,
  file = sprintf("02_RDS/mrt-yll-%s-%s.rds",
                 head(all_yrs, 1), tail(all_yrs, 1)))

## count deaths by year
subset(
  mrt,
  AGEGRP == "ALL" &
    SEX == "MF" &
    REGIOJ == "BE" &
    MEASURE == "Deaths" &
    METRIC == "Number" &
    LEVEL == 0)

## shiny output

## .. import
mrt <-
readRDS(
  sprintf("02_RDS/mrt-yll-%s-%s.rds",
          head(all_yrs, 1), tail(all_yrs, 1)))

## .. relevel SEX
mrt$SEX <-
  factor(mrt$SEX,
         levels = c("MF", "M", "F"),
         labels = c("Both sexes", "Men", "Women"))

## .. relevel REGIOJ
regioj_reg <-
  c("Brussels", "Flanders", "Wallonia", "German-speaking community")
regioj_prov <-
  c("Antwerpen", "Brabant wallon", "Hainaut", "Liege",
    "Limburg", "Luxembourg", "Namur", "Oost-Vlaanderen",
    "Vlaams-Brabant", "West-Vlaanderen")

mrt$REGIOJ <-
  factor(mrt$REGIOJ,
         levels = c("BE", "BR", "FL", "WA", "GSC",
                    "ANT", "BWA", "HAI", "LIE", "LIM",
                    "LUX", "NAM", "OVL", "VLB", "WVL"),
         labels = c("Belgium", regioj_reg, regioj_prov))

## .. add level-cause column
mrt$LEVEL_CAUSE <- paste0(mrt$LEVEL, "|", mrt$CAUSE)

## .. save as parquet
write_parquet(
  mrt,
  sink = sprintf("02_RDS/mrt-yll-%s-%s.parquet",
                 head(all_yrs, 1), tail(all_yrs, 1)))
