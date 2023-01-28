### BEBOD / SEYLL RESULTS

## extractr functions
extract <-
function(x, f) {
  y <- as.data.frame(xtabs(f, x[[1]]))
  n <- ncol(y)

  mx <- sapply(x, function(x) as.data.frame(xtabs(f, x))[, n])
  out <- cbind(y[, -n], data.frame(t(apply(mx, 1, mean_ci))))
  colnames(out) <- c(colnames(y)[-n], "mean", "2.5%", "97.5%")
  out
}

extract_std <-
function(x, f, yr, ref_be, ref_eu) {
  y <- as.data.frame(xtabs(f, x[[1]]))
  n <- ncol(y)
  yy <- y
  yy$Freq <- NULL
  yy$AGEGRP <- NULL
  yy <- yy[!duplicated(yy), , drop = FALSE]
  if (length(yy) == 0) yy <- data.frame(SEX = "MF", REGIOJ = "BE")
  if (is.null(yy$SEX)) yy$SEX <- "MF"
  if (is.null(yy$REGIOJ) & is.null(yy$PROV) & is.null(yy$GSC))
    yy$REGIOJ <- "BE"
  names(yy)[names(yy) == "PROV"] <- "REGIOJ"
  if ("GSC" %in% names(yy)) {
    yy <- subset(yy, GSC != 0)
    yy$REGIOJ <- "GSC"
    yy$GSC <- NULL
  }

  mx <- sapply(x, function(x) as.data.frame(xtabs(f, x))[, n])
  mx <- cbind(y[, -n, drop = FALSE], data.frame(mx))
  if (is.null(mx$SEX)) mx$SEX <- "MF"
  if (is.null(mx$REGIOJ) & is.null(mx$PROV) & is.null(mx$GSC))
    mx$REGIOJ <- "BE"
  names(mx)[names(mx) == "PROV"] <- "REGIOJ"
  if ("GSC" %in% names(mx)) {
    mx <- tail(mx, nrow(mx) / 2)
    mx$REGIOJ <- "GSC"
    mx$GSC <- NULL
  }
  mx <- merge(mx, subset(pop2, YEAR == yr)) ## select age/sex

  asr_be <-
  sapply(names(mx)[grepl("^X", names(mx))],
    function (j) {
      by(mx, as.list(mx[names(yy)]), 
         FUN = function(x) ASR(x[[j]], x$POP, ref_be$POP),
         simplify = FALSE)
    })
  asr_be <-
  matrix(unlist(asr_be),
         ncol = ifelse(is.null(ncol(asr_be)), length(asr_be), ncol(asr_be)))

  asr_eu <-
  sapply(names(mx)[grepl("^X", names(mx))],
    function (j) {
      by(mx, as.list(mx[names(yy)]), 
         FUN = function(x) ASR(x[[j]], x$POP, ref_eu$POP),
         simplify = FALSE)
    })
  asr_eu <-
  matrix(unlist(asr_eu),
         ncol = ifelse(is.null(ncol(asr_eu)), length(asr_be), ncol(asr_eu)))

  out_be <- cbind(yy, data.frame(t(apply(asr_be, 1, bd::mean_ci))))
  colnames(out_be) <-
    c(colnames(yy), "VAL_MEAN", "VAL_LWR", "VAL_UPR")
  out_be$AGEGRP <- "BSP"
  out_be$METRIC <- "Rate"
  out_be

  out_eu <- cbind(yy, data.frame(t(apply(asr_eu, 1, bd::mean_ci))))
  colnames(out_eu) <-
    c(colnames(yy), "VAL_MEAN", "VAL_LWR", "VAL_UPR")
  out_eu$AGEGRP <- "ESP"
  out_eu$METRIC <- "Rate"
  out_eu

  return(rbind(out_be, out_eu))
}

## insertr function
insert <-
function(x, name, fun) {
  for (i in seq_along(x)) {
    x[[i]][[name]] <- fun(x[[i]])
  }
  x
}

## export counts // all causes
export_all <-
function(x, grp) {
  ## define detailed age groups
  x <- insert(x, "AGEGRP", function(x) cut(x$AGE, grp, right = FALSE))

  ### MORTALITY

  ## extract deaths by cause, age, sex, region
  df1r <- extract(x, ~AGEGRP+SEX+REGIOJ)  # all dimensions
  df1p <- extract(x, ~AGEGRP+SEX+PROV)    # all dimensions
  df1g <- extract(x, ~AGEGRP+SEX+GSC)     # all dimensions
  df2  <- extract(x, ~AGEGRP+SEX)         # Belgium, by age & sex
  df3r <- extract(x, ~AGEGRP+REGIOJ)      # both sexes, by age & region
  df3p <- extract(x, ~AGEGRP+PROV)        # both sexes, by age & region
  df3g <- extract(x, ~AGEGRP+GSC)         # both sexes, by age & region
  df4r <- extract(x, ~SEX+REGIOJ)         # all ages, by sex & region
  df4p <- extract(x, ~SEX+PROV)           # all ages, by sex & region
  df4g <- extract(x, ~SEX+GSC)            # all ages, by sex & region
  df5  <- extract(x, ~AGEGRP)             # Belgium & both sexes, by age
  df6  <- extract(x, ~SEX)                # Belgium & all ages, by sex
  df7r <- extract(x, ~REGIOJ)             # all ages & both sexes, by reg
  df7p <- extract(x, ~PROV)               # all ages & both sexes, by reg
  df7g <- extract(x, ~GSC)                # all ages & both sexes, by reg
  df8 <- as.data.frame(t(mean_ci(sapply(x, nrow)))) # BE, ALL, MF

  ## clean GSC extracts
  df1g <- subset(df1g, GSC == 1)
  df1g$GSC <- "GSC"; names(df1g)[names(df1g) == "GSC"] <- "REGIOJ"
  df3g <- subset(df3g, GSC == 1)
  df3g$GSC <- "GSC"; names(df3g)[names(df3g) == "GSC"] <- "REGIOJ"
  df4g <- subset(df4g, GSC == 1)
  df4g$GSC <- "GSC"; names(df4g)[names(df4g) == "GSC"] <- "REGIOJ"
  df7g <- subset(df7g, GSC == 1)
  df7g$GSC <- "GSC"; names(df7g)[names(df7g) == "GSC"] <- "REGIOJ"

  ## harmonize column names to allow merge
  names(df1p)[names(df1p) == "PROV"] <- "REGIOJ"
  names(df3p)[names(df3p) == "PROV"] <- "REGIOJ"
  names(df4p)[names(df4p) == "PROV"] <- "REGIOJ"
  names(df7p)[names(df7p) == "PROV"] <- "REGIOJ"

  ## complete subsets
  df2  <- cbind(df2, REGIOJ = "BE")
  df3r <- cbind(df3r, SEX = "MF")
  df3p <- cbind(df3p, SEX = "MF")
  df3g <- cbind(df3g, SEX = "MF")
  df4r <- cbind(df4r, AGEGRP = "ALL")
  df4p <- cbind(df4p, AGEGRP = "ALL")
  df4g <- cbind(df4g, AGEGRP = "ALL")
  df5  <- cbind(df5, REGIOJ = "BE", SEX = "MF")
  df6  <- cbind(df6, REGIOJ = "BE", AGEGRP = "ALL")
  df7r <- cbind(df7r, SEX = "MF", AGEGRP = "ALL")
  df7p <- cbind(df7p, SEX = "MF", AGEGRP = "ALL")
  df7g <- cbind(df7g, SEX = "MF", AGEGRP = "ALL")
  df8  <- cbind(df8, REGIOJ = "BE", SEX = "MF", AGEGRP = "ALL")

  ## merge subsets
  df <-
    rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4r, df4p, df4g,
          df5, df6, df7r, df7p, df7g,df8)
  #str(df)

  ## reorganize columns
  df <- df[, c("AGEGRP", "SEX", "REGIOJ", "mean", "2.5%", "97.5%")]
  df_mrt <- cbind(df, MEASURE = "Deaths")


  ### YLL

  ## extract YLL by cause, age, sex, region
  df1r <- extract(x, YLL~AGEGRP+SEX+REGIOJ)  # all dimensions
  df1p <- extract(x, YLL~AGEGRP+SEX+PROV)    # all dimensions
  df1g <- extract(x, YLL~AGEGRP+SEX+GSC)     # all dimensions
  df2  <- extract(x, YLL~AGEGRP+SEX)         # Belgium, by age & sex
  df3r <- extract(x, YLL~AGEGRP+REGIOJ)      # both sexes, by age & region
  df3p <- extract(x, YLL~AGEGRP+PROV)        # both sexes, by age & region
  df3g <- extract(x, YLL~AGEGRP+GSC)         # both sexes, by age & region
  df4r <- extract(x, YLL~SEX+REGIOJ)         # all ages, by sex & region
  df4p <- extract(x, YLL~SEX+PROV)           # all ages, by sex & region
  df4g <- extract(x, YLL~SEX+GSC)            # all ages, by sex & region
  df5  <- extract(x, YLL~AGEGRP)             # Belgium & both sexes, by age
  df6  <- extract(x, YLL~SEX)                # Belgium & all ages, by sex
  df7r <- extract(x, YLL~REGIOJ)             # all ages & both sexes, by reg
  df7p <- extract(x, YLL~PROV)               # all ages & both sexes, by reg
  df7g <- extract(x, YLL~GSC)                # all ages & both sexes, by reg
  df8 <-
    as.data.frame(
      t(mean_ci(sapply(x, function(y) sum(y$YLL))))) # BE, ALL, MF

  ## clean GSC extracts
  df1g <- subset(df1g, GSC == 1)
  df1g$GSC <- "GSC"; names(df1g)[names(df1g) == "GSC"] <- "REGIOJ"
  df3g <- subset(df3g, GSC == 1)
  df3g$GSC <- "GSC"; names(df3g)[names(df3g) == "GSC"] <- "REGIOJ"
  df4g <- subset(df4g, GSC == 1)
  df4g$GSC <- "GSC"; names(df4g)[names(df4g) == "GSC"] <- "REGIOJ"
  df7g <- subset(df7g, GSC == 1)
  df7g$GSC <- "GSC"; names(df7g)[names(df7g) == "GSC"] <- "REGIOJ"

  ## harmonize column names to allow merge
  names(df1p)[names(df1p) == "PROV"] <- "REGIOJ"
  names(df3p)[names(df3p) == "PROV"] <- "REGIOJ"
  names(df4p)[names(df4p) == "PROV"] <- "REGIOJ"
  names(df7p)[names(df7p) == "PROV"] <- "REGIOJ"

  ## complete subsets
  df2  <- cbind(df2, REGIOJ = "BE")
  df3r <- cbind(df3r, SEX = "MF")
  df3p <- cbind(df3p, SEX = "MF")
  df3g <- cbind(df3g, SEX = "MF")
  df4r <- cbind(df4r, AGEGRP = "ALL")
  df4p <- cbind(df4p, AGEGRP = "ALL")
  df4g <- cbind(df4g, AGEGRP = "ALL")
  df5  <- cbind(df5, REGIOJ = "BE", SEX = "MF")
  df6  <- cbind(df6, REGIOJ = "BE", AGEGRP = "ALL")
  df7r <- cbind(df7r, SEX = "MF", AGEGRP = "ALL")
  df7p <- cbind(df7p, SEX = "MF", AGEGRP = "ALL")
  df7g <- cbind(df7g, SEX = "MF", AGEGRP = "ALL")
  df8  <- cbind(df8, REGIOJ = "BE", SEX = "MF", AGEGRP = "ALL")

  ## merge subsets
  df <-
    rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4r, df4p, df4g,
          df5, df6, df7r, df7p, df7g,df8)
  #str(df)

  ## reorganize columns
  df <- df[, c("AGEGRP", "SEX", "REGIOJ", "mean", "2.5%", "97.5%")]
  df_yll <- cbind(df, MEASURE = "YLL")


  ## export results
  rbind(df_mrt, df_yll)
}

## export age-std rates // all causes
export_std_all <-
function(x, grp, yr, ref_be, ref_eu) {
  ## define detailed age groups
  x <- insert(x, "AGEGRP", function(x) cut(x$AGE, grp, right = FALSE))

  ### MORTALITY

  ## calculate age-std death rates by sex, region
  df1r <- extract_std(x, ~AGEGRP+SEX+REGIOJ, yr, ref_be, ref_eu)
  df1p <- extract_std(x, ~AGEGRP+SEX+PROV, yr, ref_be, ref_eu)
  df1g <- extract_std(x, ~AGEGRP+SEX+GSC, yr, ref_be, ref_eu)
  df2  <- extract_std(x, ~AGEGRP+SEX, yr, ref_be, ref_eu)
  df3r <- extract_std(x, ~AGEGRP+REGIOJ, yr, ref_be, ref_eu)
  df3p <- extract_std(x, ~AGEGRP+PROV, yr, ref_be, ref_eu)
  df3g <- extract_std(x, ~AGEGRP+GSC, yr, ref_be, ref_eu)
  df4  <- extract_std(x, ~AGEGRP, yr, ref_be, ref_eu)

  ## merge subsets
  df <- rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4)
  #str(df)

  ## finalize df
  df_mrt <- cbind(df, MEASURE = "Deaths")

  ### YLL

  ## extract YLL by cause, age, sex, region
  df1r <- extract_std(x, YLL~AGEGRP+SEX+REGIOJ, yr, ref_be, ref_eu)
  df1p <- extract_std(x, YLL~AGEGRP+SEX+PROV, yr, ref_be, ref_eu)
  df1g <- extract_std(x, YLL~AGEGRP+SEX+GSC, yr, ref_be, ref_eu)
  df2  <- extract_std(x, YLL~AGEGRP+SEX, yr, ref_be, ref_eu)
  df3r <- extract_std(x, YLL~AGEGRP+REGIOJ, yr, ref_be, ref_eu)
  df3p <- extract_std(x, YLL~AGEGRP+PROV, yr, ref_be, ref_eu)
  df3g <- extract_std(x, YLL~AGEGRP+GSC, yr, ref_be, ref_eu)
  df4  <- extract_std(x, YLL~AGEGRP, yr, ref_be, ref_eu)

  ## merge subsets
  df <- rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4)
  #str(df)

  ## finalize df
  df_yll <- cbind(df, MEASURE = "YLL")

  ## export results
  rbind(df_mrt, df_yll)
}

## export counts // by cause
export <-
function(x, level, grp) {
  ## define detailed age groups
  x <- insert(x, "AGEGRP", function(x) cut(x$AGE, grp, right = FALSE))

  ## define generic name for 'cause' variable
  for (i in seq_along(x)) {
    x[[i]]$CAUSE <- x[[i]][[sprintf("GBD%s_RED4", level)]]
  }

  ### MORTALITY

  ## extract deaths by cause, age, sex, region
  df1r <- extract(x, ~CAUSE+AGEGRP+SEX+REGIOJ) # all dimensions, region
  df1p <- extract(x, ~CAUSE+AGEGRP+SEX+PROV)   # all dimensions, province
  df1g <- extract(x, ~CAUSE+AGEGRP+SEX+GSC)    # all dimensions, province
  df2  <- extract(x, ~CAUSE+AGEGRP+SEX)        # Belgium, by age & sex
  df3r <- extract(x, ~CAUSE+AGEGRP+REGIOJ)     # both sexes, by age & region
  df3p <- extract(x, ~CAUSE+AGEGRP+PROV)       # both sexes, by age & prov
  df3g <- extract(x, ~CAUSE+AGEGRP+GSC)        # both sexes, by age & prov
  df4r <- extract(x, ~CAUSE+SEX+REGIOJ)        # all ages, by sex & region
  df4p <- extract(x, ~CAUSE+SEX+PROV)          # all ages, by sex & prov
  df4g <- extract(x, ~CAUSE+SEX+GSC)           # all ages, by sex & prov
  df5  <- extract(x, ~CAUSE+AGEGRP)            # Belgium & both sexes, by age
  df6  <- extract(x, ~CAUSE+SEX)               # Belgium & all ages, by sex
  df7r <- extract(x, ~CAUSE+REGIOJ)            # all ages & sexes, by reg
  df7p <- extract(x, ~CAUSE+PROV)              # all ages & sexes, by prov
  df7g <- extract(x, ~CAUSE+GSC)               # all ages & sexes, by prov
  df8  <- extract(x, ~CAUSE)                   # Belgium, all ages & sexes

  ## clean GSC extracts
  df1g <- subset(df1g, GSC == 1)
  df1g$GSC <- "GSC"; names(df1g)[names(df1g) == "GSC"] <- "REGIOJ"
  df3g <- subset(df3g, GSC == 1)
  df3g$GSC <- "GSC"; names(df3g)[names(df3g) == "GSC"] <- "REGIOJ"
  df4g <- subset(df4g, GSC == 1)
  df4g$GSC <- "GSC"; names(df4g)[names(df4g) == "GSC"] <- "REGIOJ"
  df7g <- subset(df7g, GSC == 1)
  df7g$GSC <- "GSC"; names(df7g)[names(df7g) == "GSC"] <- "REGIOJ"

  ## harmonize column names to allow merge
  names(df1p)[names(df1p) == "PROV"] <- "REGIOJ"
  names(df3p)[names(df3p) == "PROV"] <- "REGIOJ"
  names(df4p)[names(df4p) == "PROV"] <- "REGIOJ"
  names(df7p)[names(df7p) == "PROV"] <- "REGIOJ"

  ## complete subsets
  df2  <- cbind(df2, REGIOJ = "BE")
  df3r <- cbind(df3r, SEX = "MF")
  df3p <- cbind(df3p, SEX = "MF")
  df3g <- cbind(df3g, SEX = "MF")
  df4r <- cbind(df4r, AGEGRP = "ALL")
  df4p <- cbind(df4p, AGEGRP = "ALL")
  df4g <- cbind(df4g, AGEGRP = "ALL")
  df5  <- cbind(df5, REGIOJ = "BE", SEX = "MF")
  df6  <- cbind(df6, REGIOJ = "BE", AGEGRP = "ALL")
  df7r <- cbind(df7r, SEX = "MF", AGEGRP = "ALL")
  df7p <- cbind(df7p, SEX = "MF", AGEGRP = "ALL")
  df7g <- cbind(df7g, SEX = "MF", AGEGRP = "ALL")
  df8  <- cbind(df8, REGIOJ = "BE", SEX = "MF", AGEGRP = "ALL")

  ## merge subsets
  df <-
    rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4r, df4p, df4g,
          df5, df6, df7r, df7p, df7g, df8)
  #str(df)

  ## check and remove 'Garbage code'
  if (sum(subset(df, CAUSE == "Garbage code")$mean) != 0) {
    stop("Residual garbage codes detected in ", sQuote("df_mrt"))

  } else {
    df <- subset(df, CAUSE != "Garbage code")
  }

  ## check and remove 'Stillbirth'
  if (sum(subset(df, CAUSE == "Stillbirth")$mean) != 0) {
    stop("Residual stillbirths detected in ", sQuote("df_mrt"))

  } else {
    df <- subset(df, CAUSE != "Stillbirth")
  }

  ## reorganize columns
  df <- df[, c("CAUSE", "AGEGRP", "SEX", "REGIOJ",
               "mean", "2.5%", "97.5%")]
  names(df)[names(df) == "CAUSE"] <- sprintf("CAUSE%s", level)
  df_mrt <- cbind(df, MEASURE = "Deaths")

  ### YLL

  ## extract YLL by cause, age, sex, region
  df1r <- extract(x, YLL~CAUSE+AGEGRP+SEX+REGIOJ) # all dimensions, by region
  df1p <- extract(x, YLL~CAUSE+AGEGRP+SEX+PROV)   # all dimensions, by prov
  df1g <- extract(x, YLL~CAUSE+AGEGRP+SEX+GSC)    # all dimensions, by prov
  df2  <- extract(x, YLL~CAUSE+AGEGRP+SEX)        # BE, by age & sex
  df3r <- extract(x, YLL~CAUSE+AGEGRP+REGIOJ)     # both sex, by age & region
  df3p <- extract(x, YLL~CAUSE+AGEGRP+PROV)       # both sex, by age & prov
  df3g <- extract(x, YLL~CAUSE+AGEGRP+GSC)        # both sex, by age & prov
  df4r <- extract(x, YLL~CAUSE+SEX+REGIOJ)        # all age, by sex & region
  df4p <- extract(x, YLL~CAUSE+SEX+PROV)          # all age, by sex & prov
  df4g <- extract(x, YLL~CAUSE+SEX+GSC)           # all age, by sex & prov
  df5  <- extract(x, YLL~CAUSE+AGEGRP)            # BE & both sexes, by age
  df6  <- extract(x, YLL~CAUSE+SEX)               # BE & all age, by sex
  df7r <- extract(x, YLL~CAUSE+REGIOJ)            # all age & sexes, by reg
  df7p <- extract(x, YLL~CAUSE+PROV)              # all age & sexes, by prov
  df7g <- extract(x, YLL~CAUSE+GSC)               # all age & sexes, by prov
  df8  <- extract(x, YLL~CAUSE)                   # BE, all age & sexes

  ## clean GSC extracts
  df1g <- subset(df1g, GSC == 1)
  df1g$GSC <- "GSC"; names(df1g)[names(df1g) == "GSC"] <- "REGIOJ"
  df3g <- subset(df3g, GSC == 1)
  df3g$GSC <- "GSC"; names(df3g)[names(df3g) == "GSC"] <- "REGIOJ"
  df4g <- subset(df4g, GSC == 1)
  df4g$GSC <- "GSC"; names(df4g)[names(df4g) == "GSC"] <- "REGIOJ"
  df7g <- subset(df7g, GSC == 1)
  df7g$GSC <- "GSC"; names(df7g)[names(df7g) == "GSC"] <- "REGIOJ"

  ## harmonize column names to allow merge
  names(df1p)[names(df1p) == "PROV"] <- "REGIOJ"
  names(df3p)[names(df3p) == "PROV"] <- "REGIOJ"
  names(df4p)[names(df4p) == "PROV"] <- "REGIOJ"
  names(df7p)[names(df7p) == "PROV"] <- "REGIOJ"

  ## complete subsets
  df2  <- cbind(df2, REGIOJ = "BE")
  df3r <- cbind(df3r, SEX = "MF")
  df3p <- cbind(df3p, SEX = "MF")
  df3g <- cbind(df3g, SEX = "MF")
  df4r <- cbind(df4r, AGEGRP = "ALL")
  df4p <- cbind(df4p, AGEGRP = "ALL")
  df4g <- cbind(df4g, AGEGRP = "ALL")
  df5  <- cbind(df5, REGIOJ = "BE", SEX = "MF")
  df6  <- cbind(df6, REGIOJ = "BE", AGEGRP = "ALL")
  df7r <- cbind(df7r, SEX = "MF", AGEGRP = "ALL")
  df7p <- cbind(df7p, SEX = "MF", AGEGRP = "ALL")
  df7g <- cbind(df7g, SEX = "MF", AGEGRP = "ALL")
  df8  <- cbind(df8, REGIOJ = "BE", SEX = "MF", AGEGRP = "ALL")

  ## merge subsets
  df <-
    rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4r, df4p, df4g,
          df5, df6, df7r, df7p, df7g, df8)
  #str(df)

  ## check and remove 'Garbage code'
  if (sum(subset(df, CAUSE == "Garbage code")$mean) != 0) {
    stop("Residual garbage codes detected in ", sQuote("df_yll"))

  } else {
    df <- subset(df, CAUSE != "Garbage code")
  }

  ## check and remove 'Stillbirth'
  if (sum(subset(df, CAUSE == "Stillbirth")$mean) != 0) {
    stop("Residual stillbirths detected in ", sQuote("df_yll"))

  } else {
    df <- subset(df, CAUSE != "Stillbirth")
  }

  ## reorganize columns
  df <- df[, c("CAUSE", "AGEGRP", "SEX", "REGIOJ", "mean", "2.5%", "97.5%")]
  names(df)[names(df) == "CAUSE"] <- sprintf("CAUSE%s", level)
  df_yll <- cbind(df, MEASURE = "YLL")


  ## combine and export results
  rbind(df_mrt, df_yll)
}


## export age-std rates // by cause
export_std <-
function(x, level, grp, yr, ref_be, ref_eu) {
  ## define detailed age groups
  x <- insert(x, "AGEGRP", function(x) cut(x$AGE, grp, right = FALSE))

  ## define generic name for 'cause' variable
  for (i in seq_along(x)) {
    x[[i]]$CAUSE <- x[[i]][[sprintf("GBD%s_RED4", level)]]
  }

  ### MORTALITY

  ## extract deaths by cause, age, sex, region
  df1r <- extract_std(x, ~CAUSE+AGEGRP+SEX+REGIOJ, yr, ref_be, ref_eu)
  df1p <- extract_std(x, ~CAUSE+AGEGRP+SEX+PROV, yr, ref_be, ref_eu)
  df1g <- extract_std(x, ~CAUSE+AGEGRP+SEX+GSC, yr, ref_be, ref_eu)
  df2  <- extract_std(x, ~CAUSE+AGEGRP+SEX, yr, ref_be, ref_eu)
  df3r <- extract_std(x, ~CAUSE+AGEGRP+REGIOJ, yr, ref_be, ref_eu)
  df3p <- extract_std(x, ~CAUSE+AGEGRP+PROV, yr, ref_be, ref_eu)
  df3g <- extract_std(x, ~CAUSE+AGEGRP+GSC, yr, ref_be, ref_eu)
  df4  <- extract_std(x, ~CAUSE+AGEGRP, yr, ref_be, ref_eu)

  ## merge subsets
  df <- rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4)
  #str(df)

  ## check and remove 'Garbage code'
  if (sum(subset(df, CAUSE == "Garbage code")$mean) != 0) {
    stop("Residual garbage codes detected in ", sQuote("df_mrt"))

  } else {
    df <- subset(df, CAUSE != "Garbage code")
  }

  ## check and remove 'Stillbirth'
  if (sum(subset(df, CAUSE == "Stillbirth")$mean) != 0) {
    stop("Residual stillbirths detected in ", sQuote("df_mrt"))

  } else {
    df <- subset(df, CAUSE != "Stillbirth")
  }

  ## reorganize columns
  names(df)[names(df) == "CAUSE"] <- sprintf("CAUSE%s", level)
  df_mrt <- cbind(df, MEASURE = "Deaths")

  ### YLL

  ## extract YLL by cause, age, sex, region
  df1r <- extract_std(x, YLL~CAUSE+AGEGRP+SEX+REGIOJ, yr, ref_be, ref_eu)
  df1p <- extract_std(x, YLL~CAUSE+AGEGRP+SEX+PROV, yr, ref_be, ref_eu)
  df1g <- extract_std(x, YLL~CAUSE+AGEGRP+SEX+GSC, yr, ref_be, ref_eu)
  df2  <- extract_std(x, YLL~CAUSE+AGEGRP+SEX, yr, ref_be, ref_eu)
  df3r <- extract_std(x, YLL~CAUSE+AGEGRP+REGIOJ, yr, ref_be, ref_eu)
  df3p <- extract_std(x, YLL~CAUSE+AGEGRP+PROV, yr, ref_be, ref_eu)
  df3g <- extract_std(x, YLL~CAUSE+AGEGRP+GSC, yr, ref_be, ref_eu)
  df4  <- extract_std(x, YLL~CAUSE+AGEGRP, yr, ref_be, ref_eu)

  ## merge subsets
  df <- rbind(df1r, df1p, df1g, df2, df3r, df3p, df3g, df4)
  #str(df)

  ## check and remove 'Garbage code'
  if (sum(subset(df, CAUSE == "Garbage code")$mean) != 0) {
    stop("Residual garbage codes detected in ", sQuote("df_yll"))

  } else {
    df <- subset(df, CAUSE != "Garbage code")
  }

  ## check and remove 'Stillbirth'
  if (sum(subset(df, CAUSE == "Stillbirth")$mean) != 0) {
    stop("Residual stillbirths detected in ", sQuote("df_yll"))

  } else {
    df <- subset(df, CAUSE != "Stillbirth")
  }

  ## reorganize columns
  names(df)[names(df) == "CAUSE"] <- sprintf("CAUSE%s", level)
  df_yll <- cbind(df, MEASURE = "YLL")

  ## combine and export results
  rbind(df_mrt, df_yll)
}

summarize_all <-
function(mrt_sim, level, yr, pop, pop2, grp, grp2) {
  ## get counts for all combinations
  mrt_out1 <-
  switch(
    as.character(level),
    "0" = export_all(mrt_sim, grp),
    "1" = export(mrt_sim, level, grp),
    "2" = export(mrt_sim, level, grp),
    "3" = export(mrt_sim, level, grp))
  mrt_out1$METRIC <- "Number"

  ## change colnames
  names(mrt_out1)[names(mrt_out1) == "mean"]  <- "VAL_MEAN"
  names(mrt_out1)[names(mrt_out1) == "2.5%"]  <- "VAL_LWR"
  names(mrt_out1)[names(mrt_out1) == "97.5%"] <- "VAL_UPR"

  ## add rates
  mrt_out1_rt <- merge(mrt_out1, subset(pop, YEAR == yr))
  mrt_out1_rt$VAL_MEAN  <- 1e5 * mrt_out1_rt$VAL_MEAN / mrt_out1_rt$POP
  mrt_out1_rt$VAL_LWR <- 1e5 * mrt_out1_rt$VAL_LWR  / mrt_out1_rt$POP
  mrt_out1_rt$VAL_UPR <- 1e5 * mrt_out1_rt$VAL_UPR / mrt_out1_rt$POP
  mrt_out1_rt$METRIC <- "Rate"

  ## merge numbers and rates
  mrt_out1 <- rbind(mrt_out1, mrt_out1_rt[, colnames(mrt_out1)])

  ## get age-std rates for all combinations
  mrt_out2 <-
  switch(
    as.character(level),
    "0" = export_std_all(mrt_sim, grp2, yr, ref_be, ref_eu),
    "1" = export_std(mrt_sim, level, grp2, yr, ref_be, ref_eu),
    "2" = export_std(mrt_sim, level, grp2, yr, ref_be, ref_eu),
    "3" = export_std(mrt_sim, level, grp2, yr, ref_be, ref_eu))

  ## merge datasets
  mrt_out <- merge(mrt_out1, mrt_out2, all = TRUE)
  #str(mrt_out)

  ## add year
  mrt_out$YEAR <- yr

  ## add cause hierarchy
  if (level == 3) {  # add level 1 and 2
    mrt_out$CAUSE2 <-
      BeBOD:::causelist$Level2[
        match(mrt_out$CAUSE3, BeBOD:::causelist$Level3)]
    mrt_out$CAUSE1 <-
      BeBOD:::causelist$Level1[
        match(mrt_out$CAUSE3, BeBOD:::causelist$Level3)]
  }
  if (level == 2) {  # add level 1
    mrt_out$CAUSE1 <-
      BeBOD:::causelist$Level1[
        match(mrt_out$CAUSE2, BeBOD:::causelist$Level2)]
  }

  ## reorder columns
  cols <-
  switch(
    as.character(level),
    "0" = NULL,
    "1" = "CAUSE1",
    "2" = c("CAUSE2", "CAUSE1"),
    "3" = c("CAUSE3", "CAUSE2", "CAUSE1"))

  mrt_out <-
    mrt_out[, c(cols,
                "YEAR", "AGEGRP", "SEX", "REGIOJ", "MEASURE", "METRIC",
                "VAL_MEAN", "VAL_LWR", "VAL_UPR")]

  ## write output
  saveRDS(
    mrt_out,
    file = sprintf("02_RDS/mrt-yll-level%s-%s.rds", level, yr))
}

## population data // export
#pop  <- read.csv2("../../../03_POP/BE-POP-2004-2019-v2.csv")
#names(pop)[names(pop) == "AGE"] <- "AGEGRP"
#names(pop)[names(pop) == "REGION"] <- "REGIOJ"

## population data // BE standard
pop2 <- read.csv2("../../../03_POP/BE-POP-2004-2020-5y-prov.csv")
names(pop2)[names(pop2) == "AGE"] <- "AGEGRP"
names(pop2)[names(pop2) == "REGION"] <- "REGIOJ"
pop2$AGEGRP2 <-
suppressWarnings(
sprintf("[%s,%s)",
        as.numeric(gsub("X", "", pop2$AGEGRP)),
        as.numeric(gsub("X", "", pop2$AGEGRP)) + 5))
pop2$AGEGRP2[pop2$AGEGRP == "ALL"] <- "ALL"
pop2$AGEGRP2[pop2$AGEGRP2 == "[85,90)"] <- "[85,Inf)"
pop2$AGEGRP <- pop2$AGEGRP2
pop2$AGEGRP2 <- NULL

pop <- pop2

## define BE population reference (most recent year, total country)
(ref_yr <- max(all_yrs))
ref_be <-
subset(pop2,
  REGIOJ == "BE" & SEX == "MF" & AGEGRP != "ALL" & YEAR == ref_yr)
ref_be <- ref_be[, c("AGEGRP", "POP")]

## define EU population reference (ESP2013)
ref_eu <- read.csv("../../../03_POP/EU_POP.csv")
names(ref_eu) <- c("AGEGRP", "POP")
ref_eu$AGEGRP <-
sprintf("[%s,%s)",
        as.numeric(gsub("X", "", ref_eu$AGEGRP)),
        as.numeric(gsub("X", "", ref_eu$AGEGRP)) + 5)
ref_eu[ref_eu == "[85,90)"] <- "[85,Inf)"

## calculate age-std rates
ASR <-
function(nr, pop, pop_ref) {
  rt_age <- nr / pop
  rt_std <- 1e5 * weighted.mean(rt_age, pop_ref)
  rt_std
}

## summarize results for different years

# .. define age groups
grp  <- #c(0, 5, 15, 45, 65, 85, Inf)  # for export
grp2 <- c(seq(0, 85, 5), Inf)         # for age-std rates

# .. run all years & levels
options(future.globals.maxSize = Inf) # disable check
plan(multisession(workers = 8))
future_sapply(all_yrs,
  function(y) {
    mrt_sim_y <- readRDS(sprintf("01_SIM/mrt-sim-%s.rds", y))
    for (l in 0:3)
      summarize_all(
        mrt_sim = mrt_sim_y, level = l, yr = y, pop, pop2, grp, grp2)
    })
