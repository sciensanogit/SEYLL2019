#' ---
#' title: BEBOD SEYLL 2019
#' subtitle: create dataset
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#' ---

#' # Data

## read data
mrt_all <- 
readRDS(paste0("//sciensano.be/FS/1140_DATA/Causes_of_death/Users/BD/",
               "dth_2000_2019_users.rds"))
str(mrt_all)


#' # Clean data

## check for stillbirths .. should be 0 .. if not replace by NA
sum(grepl("P95", mrt_all$ICD))
sum(grepl("P95", mrt_all$UCAUSE4))

sum(grepl("P95", mrt_all$CD_MLTPL_COD))
subset(mrt_all, grepl("P95", mrt_all$CD_MLTPL_COD))
mrt_all[mrt_all$CLE == "0108A195", "CD_MLTPL_COD"] <-
  "P017 P026 P034 P209"

sum(grepl("P95", mrt_all$MIIICVC4A))
subset(mrt_all, grepl("P95", mrt_all$MIIICVC4A))[, "MIIICVC4A"]
mrt_all[grepl("P95", mrt_all$MIIICVC4A), "MIIICVC4A"] <- NA
sum(grepl("P95", mrt_all$MIIICVC4B))  # should be 0
sum(grepl("P95", mrt_all$MIIICVC4C))  # should be 0
sum(grepl("P95", mrt_all$MIIICVC4D))  # should be 0
sum(grepl("P95", mrt_all$MIIICVC4E))  # should be 0
sum(grepl("P95", mrt_all$MIIICVC4F))  # should be 0
sum(grepl("P95", mrt_all$MIIICVC4G))  # should be 0

## check for SIDS in non-newborns .. should be age <= 2
table(subset(mrt_all, grepl("R95", ICD))$AGE)
table(subset(mrt_all, grepl("R95", UCAUSE4))$AGE)
table(subset(mrt_all, grepl("R95", MIIICVC4A))$AGE)
table(subset(mrt_all, grepl("R95", MIIICVC4B))$AGE)
table(subset(mrt_all, grepl("R95", MIIICVC4C))$AGE)
table(subset(mrt_all, grepl("R95", MIIICVC4D))$AGE)
table(subset(mrt_all, grepl("R95", MIIICVC4E))$AGE) #*
table(subset(mrt_all, grepl("R95", MIIICVC4F))$AGE) #*
table(subset(mrt_all, grepl("R95", MIIICVC4G))$AGE) #*

subset(mrt_all, grepl("R95", MIIICVC4A) & AGE > 5)
id <- grepl("R95", mrt_all$MIIICVC4A) & mrt_all$AGE > 5; sum(id)
mrt_all[id, "MIIICVC4A"] <- "R96"

subset(mrt_all, grepl("R95", MIIICVC4D) & AGE > 5)
id <- grepl("R95", mrt_all$MIIICVC4D) & mrt_all$AGE > 5; sum(id)
mrt_all[id, "MIIICVC4D"] <- "R96"

subset(mrt_all, grepl("R95", MIIICVC4E) & AGE > 5)
id <- grepl("R95", mrt_all$MIIICVC4E) & mrt_all$AGE > 5; sum(id)
mrt_all[id, "MIIICVC4E"] <- "R96"

## check if 'ICD' not equal to 'UCAUSE4'
mrt_wrong_id <- mrt_all$ICD != substr(mrt_all$UCAUSE4, 0, 3)
mrt_wrong <- mrt_all[mrt_wrong_id, ]
mrt_wrong[, c("YEAR", "AGE", "ICD", "UCAUSE4")]

#write_cb(cbind(
#  map_gbd(mrt_wrong$ICD)[, "cause4"],
#  map_gbd(mrt_wrong$UCAUSE4)[, "cause4"]))

## .. use ICD unless it is a garbage code
mrt_wrong_ucause4 <- mrt_all[mrt_wrong_id, "ICD"]
mrt_wrong_ucause4[
  map_gbd(mrt_wrong$ICD)[, "cause4"] == "Garbage code" &
    !is.na(map_gbd(mrt_wrong$UCAUSE4)[, "cause4"]) &
    map_gbd(mrt_wrong$UCAUSE4)[, "cause4"] != "Garbage code"] <-
mrt_all[mrt_wrong_id, "UCAUSE4"][
  map_gbd(mrt_wrong$ICD)[, "cause4"] == "Garbage code" &
    !is.na(map_gbd(mrt_wrong$UCAUSE4)[, "cause4"]) &
    map_gbd(mrt_wrong$UCAUSE4)[, "cause4"] != "Garbage code"]

mrt_all[mrt_wrong_id, "UCAUSE4"] <- mrt_wrong_ucause4

## check for missing UCAUSE4
subset(mrt_all, UCAUSE4 %in% c("", "*"))  # empty?

## check for U codes (provisional codes)
## .. if possible, define as specific codes in 'ICD-GBD-mapping'
subset(mrt_all, grepl("U", UCAUSE4))
subset(mrt_all, grepl("U", MIIICVC4D))
table(subset(mrt_all, grepl("U", MIIICVC4C))$MIIICVC4C)
table(subset(mrt_all, grepl("U", MIIICVC4B))$MIIICVC4B)
table(subset(mrt_all, grepl("U", MIIICVC4A))$MIIICVC4A)
subset(mrt_all, MIIICVC4A == "U4690")

table(subset(mrt_all, grepl("U", MIIICVC4E))$MIIICVC4E) #*
table(subset(mrt_all, grepl("U", MIIICVC4F))$MIIICVC4F) #*
table(subset(mrt_all, grepl("U", MIIICVC4G))$MIIICVC4G) #*

## remove wrong entries as cause
mrt_all$MIIICVC2[mrt_all$MIIICVC2 == "*"] <- NA
mrt_all$MIIICVC2[mrt_all$MIIICVC2 == "0"] <- NA

mrt_all$MIIICVC2 <- gsub(",", "", mrt_all$MIIICVC2)
mrt_all$MIIICVC4A <- gsub(",", "", mrt_all$MIIICVC4A)
mrt_all$MIIICVC4B <- gsub(",", "", mrt_all$MIIICVC4B)
mrt_all$MIIICVC4C <- gsub(",", "", mrt_all$MIIICVC4C)
mrt_all$MIIICVC4D <- gsub(",", "", mrt_all$MIIICVC4D)
mrt_all$MIIICVC4E <- gsub(",", "", mrt_all$MIIICVC4E) #*
mrt_all$MIIICVC4F <- gsub(",", "", mrt_all$MIIICVC4F) #*
mrt_all$MIIICVC4G <- gsub(",", "", mrt_all$MIIICVC4G) #*

mrt_all$MIIICVC2 <- toupper(mrt_all$MIIICVC2)
mrt_all$MIIICVC4A <- toupper(mrt_all$MIIICVC4A)
mrt_all$MIIICVC4B <- toupper(mrt_all$MIIICVC4B)
mrt_all$MIIICVC4C <- toupper(mrt_all$MIIICVC4C)
mrt_all$MIIICVC4D <- toupper(mrt_all$MIIICVC4D)
mrt_all$MIIICVC4E <- toupper(mrt_all$MIIICVC4E) #*
mrt_all$MIIICVC4F <- toupper(mrt_all$MIIICVC4F) #*
mrt_all$MIIICVC4G <- toupper(mrt_all$MIIICVC4G) #*


#' # Descriptive

## deaths by year
table(mrt_all$YEAR)


#' # Define variables

## redefine region
mrt_all$REGIOJ[mrt_all$REGIOJ == 1] <- "FL"
mrt_all$REGIOJ[mrt_all$REGIOJ == 2] <- "BR"
mrt_all$REGIOJ[mrt_all$REGIOJ == 3] <- "WA"
table(mrt_all$REGIOJ)

## redefine province // take out 12=BCR
mrt_all$PROV <-
factor(
  mrt_all$PROV,
  c(1, 3:11),
  c("ANT", "WVL", "OVL", "HAI", "LIE", "LIM", "LUX", "NAM", "VLB", "BWA"))
table(mrt_all$PROV, useNA = "always")

## rename GSC
names(mrt_all)[names(mrt_all) == "SUBPROV_RES"] <- "GSC"

## redefine sex
mrt_all$SEX <- factor(mrt_all$SEX, levels = 1:2, labels = c("M", "F"))
table(mrt_all$SEX)

## create AGE CAT
mrt_all$AGEGRP <-
  cut(mrt_all$AGE, c(0, 5, 15, 45, 65, 85, Inf), right = FALSE)
mrt_all$AGEGRP <- factor(mrt_all$AGEGRP)

## .. take out observations with missing age or sex
range(mrt_all$AGE, na.rm = TRUE)
table(mrt_all$AGEGRP, useNA = "always")
table(mrt_all$SEX, useNA = "always")

mrt_all <- subset(mrt_all, !is.na(AGEGRP))
mrt_all <- subset(mrt_all, !is.na(SEX))

## create AGESEX
mrt_all$AGEGRPSEX <- paste(mrt_all$AGEGRP, mrt_all$SEX)
mrt_all$AGEGRPSEX <- factor(mrt_all$AGEGRPSEX, unique(mrt_all$AGEGRPSEX))
table(mrt_all$AGEGRPSEX, useNA = "always")  # should have no NA values

## calculate YLLs
mrt_all$YLL <- rsle(mrt_all$AGE)
sum(is.na(mrt_all$AGE))        # should be 0
mrt_all[is.na(mrt_all$AGE), ]  # empty?
sum(is.na(mrt_all$YLL))        # should be 0
mrt_all[is.na(mrt_all$YLL), ]  # empty?


#' ## GBD causes

## (1) map to GBD cause list
in_gbd1 <- mrt_all$UCAUSE4 %in% gbd$icd_code
table(in_gbd1)
sort(table(mrt_all$UCAUSE4[!in_gbd1]))

## drop 4th digit from non-mapped causes
mrt_all$UCAUSE4_2 <- mrt_all$UCAUSE4
mrt_all$UCAUSE4_2[!in_gbd1] <- substr(mrt_all$UCAUSE4_2[!in_gbd1], 1, 3)

## (2) map to GBD cause list
in_gbd2 <- mrt_all$UCAUSE4_2 %in% gbd$icd_code
table(in_gbd2)
sort(table(mrt_all$UCAUSE4_2[!in_gbd2]))

## define GBD cause
mrt_all$CAUSE4 <- 
  gbd$yll_cause_name[match(mrt_all$UCAUSE4_2, gbd$icd_code)]
sum(is.na(mrt_all$CAUSE4))       # zero?
mrt_all[is.na(mrt_all$CAUSE4), ] # empty?
tail(sort(table(mrt_all$CAUSE4)))

## add other cause levels

## check if all causes are in mastercauselist -> should be zero!
## .. if not, add manually to BeBOD/data-raw/mastercauselist.xlsx
id <- mrt_all$CAUSE4 %in% causelist$Level4
sort(table(mrt_all$CAUSE4[!id]))

mrt_all$CAUSE3 <-
  causelist$Level3[match(mrt_all$CAUSE4, causelist$Level4)]
tail(sort(table(mrt_all$CAUSE3)))

mrt_all$CAUSE2 <-
  causelist$Level2[match(mrt_all$CAUSE4, causelist$Level4)]
tail(sort(table(mrt_all$CAUSE2)))

mrt_all$CAUSE1 <-
  causelist$Level1[match(mrt_all$CAUSE4, causelist$Level4)]
sort(table(mrt_all$CAUSE1))

## add ICD10 names
id <- match(mrt_all$UCAUSE4_2, gbd$icd_code)
mrt_all$ICD_NAME <- gbd$icd_name[id]
sum(is.na(mrt_all$ICD_NAME))
sort(table(subset(mrt_all, is.na(ICD_NAME))$UCAUSE4_2))

## define list of multiple causes of death
## .. compile ABCD causes
MCOD_ICD_ABCD <-
mrt_all[, c("MIIICVC2", "MIIICVC4A", "MIIICVC4B", "MIIICVC4C", "MIIICVC4D",
            "MIIICVC4E", "MIIICVC4F", "MIIICVC4G")]
MCOD_ICD_ABCD <- sapply(MCOD_ICD_ABCD, substr, 0, 4)
MCOD_ICD_ABCD <- as.list(as.data.frame(t(MCOD_ICD_ABCD)))
MCOD_ICD_ABCD <- lapply(MCOD_ICD_ABCD, function(x) x[x != ""])
mrt_all$MCOD_ICD <- MCOD_ICD_ABCD

## .. clean CD_MLTPL_COD
mrt_all$CD_MLTPL_COD <- gsub("\\|", "", mrt_all$CD_MLTPL_COD)
mrt_all$CD_MLTPL_COD <- gsub("\\(", "", mrt_all$CD_MLTPL_COD)
mrt_all$CD_MLTPL_COD <- gsub("/", " ", mrt_all$CD_MLTPL_COD)
mrt_all$CD_MLTPL_COD <- gsub("\\s+", " ", trimws(mrt_all$CD_MLTPL_COD))

## .. compile into unique codes
MCOD_ICD_MLTPL <- strsplit(mrt_all$CD_MLTPL_COD, " ")
table(sapply(MCOD_ICD_MLTPL, length))
table(sapply(unlist(MCOD_ICD_MLTPL), nchar))
table(unlist(MCOD_ICD_MLTPL))
MCOD_ICD_MLTPL[sapply(MCOD_ICD_MLTPL, length) > 20]
mrt_all$MCOD_ICD[mrt_all$CD_MLTPL_COD != ""] <-
  MCOD_ICD_MLTPL[sapply(MCOD_ICD_MLTPL, length) != 0]

## .. trim codes
f4 <-
function(i, j) {
  i[j] <- substr(i[j], 1, 4)
  return(unlist(i))
}

f3 <-
function(i, j) {
  i[j] <- substr(i[j], 1, 3)
  return(unlist(i))
}

not_gbd1 <- lapply(mrt_all$MCOD_ICD, function(x) is.na(fmatch(x, gbd[[3]])))
table(unlist(not_gbd1))
sort(table(unlist(mapply(`[`, mrt_all$MCOD_ICD, not_gbd1))))

b <- mapply(f4, i = mrt_all$MCOD_ICD, j = not_gbd1)
not_gbd2 <- lapply(b, function(x) is.na(fmatch(x, gbd[[3]])))
table(unlist(not_gbd2))
sort(table(unlist(mapply(`[`, b, not_gbd2))))

c <- mapply(f3, i = b, j = not_gbd2)
not_gbd3 <- lapply(c, function(x) is.na(fmatch(x, gbd[[3]])))
table(unlist(not_gbd3))
sort(table(unlist(mapply(`[`, c, not_gbd3))))

mrt_all$MCOD_ICD_ORG <- mrt_all$MCOD_ICD
mrt_all$MCOD_ICD <- c

## .. create individual variables
max(sapply(mrt_all$MCOD_ICD, length)) # check if we need more vars!!
table(sapply(mrt_all$MCOD_ICD, length))
mrt_all$MCOD_ICD_01 <- sapply(mrt_all$MCOD_ICD, function(x) x[1])
mrt_all$MCOD_ICD_02 <- sapply(mrt_all$MCOD_ICD, function(x) x[2])
mrt_all$MCOD_ICD_03 <- sapply(mrt_all$MCOD_ICD, function(x) x[3])
mrt_all$MCOD_ICD_04 <- sapply(mrt_all$MCOD_ICD, function(x) x[4])
mrt_all$MCOD_ICD_05 <- sapply(mrt_all$MCOD_ICD, function(x) x[5])
mrt_all$MCOD_ICD_06 <- sapply(mrt_all$MCOD_ICD, function(x) x[6])
mrt_all$MCOD_ICD_07 <- sapply(mrt_all$MCOD_ICD, function(x) x[7])
mrt_all$MCOD_ICD_08 <- sapply(mrt_all$MCOD_ICD, function(x) x[8])
mrt_all$MCOD_ICD_09 <- sapply(mrt_all$MCOD_ICD, function(x) x[9])
mrt_all$MCOD_ICD_10 <- sapply(mrt_all$MCOD_ICD, function(x) x[10])
mrt_all$MCOD_ICD_11 <- sapply(mrt_all$MCOD_ICD, function(x) x[11])
mrt_all$MCOD_ICD_12 <- sapply(mrt_all$MCOD_ICD, function(x) x[12])
mrt_all$MCOD_ICD_13 <- sapply(mrt_all$MCOD_ICD, function(x) x[13])
mrt_all$MCOD_ICD_14 <- sapply(mrt_all$MCOD_ICD, function(x) x[14])
mrt_all$MCOD_ICD_15 <- sapply(mrt_all$MCOD_ICD, function(x) x[15])
mrt_all$MCOD_ICD_16 <- sapply(mrt_all$MCOD_ICD, function(x) x[16])
mrt_all$MCOD_ICD_17 <- sapply(mrt_all$MCOD_ICD, function(x) x[17])
mrt_all$MCOD_ICD_18 <- sapply(mrt_all$MCOD_ICD, function(x) x[18])
mrt_all$MCOD_ICD_19 <- sapply(mrt_all$MCOD_ICD, function(x) x[19])
mrt_all$MCOD_ICD_20 <- sapply(mrt_all$MCOD_ICD, function(x) x[20])
mrt_all$MCOD_ICD_21 <- sapply(mrt_all$MCOD_ICD, function(x) x[21])
mrt_all$MCOD_ICD_22 <- sapply(mrt_all$MCOD_ICD, function(x) x[22])
mrt_all$MCOD_ICD_23 <- sapply(mrt_all$MCOD_ICD, function(x) x[23])

## define list of *specific* multiple causes of death
MCOD_GBD <-
  lapply(mrt_all$MCOD_ICD, function(x) gbd$yll_cause[fmatch(x, gbd[[3]])])
MCOD_GBD <- sapply(MCOD_GBD, paste, collapse = "|")
MCOD_GBD <- gsub("NA", "_gc", MCOD_GBD)
MCOD_GBD <- strsplit(MCOD_GBD, "\\|")
mrt_all$MCOD_GBD <- MCOD_GBD

## check if all causes in causelist
## .. add these to BeBOD/mastercauselist.xlsx
id <- !(unlist(MCOD_GBD) %fin% causelist$Code)
unique(unlist(MCOD_GBD)[id])
unique(unlist(mrt_all$MCOD_ICD)[id])

## check if any "_u" codes
id <- sapply(MCOD_GBD, function(x) any(x %fin% "_u"))
mrt_all$CD_MTLPL_COD[id]

## check if any "_sb" codes
id <- sapply(mrt_all$MCOD_ICD, function(x) any(substr(x, 1, 3) == "P95"))
mrt_all$CD_MTLPL_COD[id]


#' ## Explore IDD

sum(mrt_all$CAUSE4 == "Garbage code")
sum(mrt_all$CAUSE4 == "Garbage code") / nrow(mrt_all)

boxplot(AGE ~ CAUSE4 == "Garbage code", mrt_all)
with(mrt_all, tapply(AGE, CAUSE4 == "Garbage code", mean))

mosaicplot(SEX ~ CAUSE4 == "Garbage code", mrt_all)
with(mrt_all, tapply(CAUSE4 == "Garbage code", SEX, mean))
with(mrt_all, tapply(SEX == "M", CAUSE4 == "Garbage code", mean))
with(mrt_all, tapply(SEX == "F", CAUSE4 == "Garbage code", mean))

mosaicplot(REGIOJ ~ CAUSE4 == "Garbage code", mrt_all)
tab_reg <-
  addmargins(with(mrt_all, table(REGIOJ, CAUSE4 == "Garbage code")),
             margin = 1)
t(t(tab_reg[1:3, ]) / tab_reg[4, ])
tab_reg <-
  addmargins(with(mrt_all, table(REGIOJ, CAUSE4 == "Garbage code")),
             margin = 2)
tab_reg[, 1:2] / tab_reg[, 3]

## show non-defined garbage codes
## .. need to add these manually to 'BeBOD-IDD'
## .. check definitions in 'intermediate GC' and 'master_cause_list'
## .. for packages, check Johnson et al 2021
all_gbd <- map_gbd(mrt_all$UCAUSE4_2)[, "cause4"]
sum(is.na(all_gbd))  # zero?

all_gbd_icd <- unique(mrt_all[all_gbd == "Garbage code", "UCAUSE4_2"])
all(all_gbd_icd %in% idd$ICD_CODE)
(add_idd <- sort(all_gbd_icd[!(all_gbd_icd %in% idd$ICD_CODE)]))
#write_cb(cbind(add_idd, gbd$icd_name[match(add_idd, gbd$icd_code)]))

# U-codes
# u <- sapply(mrt_all$MCOD_ICD, function(x) any(grepl("^U", x)))
# bd::write_cb(mrt_all[u, 1:46], col.names = TRUE)

##rmarkdown::render("01-SEYLL-create.R")