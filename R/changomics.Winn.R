# WiNN version 0.4
library(hwwntest)
library(gam)
library(mgcv)
library(stringr)
library(lawstat) # <= Form Leveve test. Temporary until we get the "car" library working


#####################################################
# Utility function for test.wn
# Find greatest power of 2 that is closest to "n"
# ("closest.pw2", could be > or < "n") and the greatest
# power of 2 that is closest to "n" and < "n"
# ("greatest.pw2")
#####################################################
greatest.closest.power.2 <- function(n) {
  i <- 0

  repeat {
    if (n >= 2 ^ i & n < 2 ^ (i + 1)) {
      break
    }
    else{
      i <- i + 1
    }
  }

  i.closest <- i

  # Is this really the power od 2 closest to "n"?
  # (n=255 would return a "power of 2" i value = 7, while it is closer to 2^8)
  if (abs(n - 2 ^ (i + 1)) < (n - 2 ^ i)) {
    i.closest <- (i + 1)
  }

  return(list(closest.pw2 = i.closest, greatest.pw2 = i))
}


######################################
# White Noise test
# The test combines hwwntest and Box
# tests
######################################
test.wn <- function(x, box.test.lag = 20) {
  # The p-value for WN to return
  p <- NA

  # Remove NAs
  x.no.na <- x[!is.na(x)]

  # The WN test can be applied for series lenghts >= 16
  if (length(x.no.na) >= 16) {
    # Determine greatest power of 2 which is < length of data
    z <- greatest.closest.power.2(length(x.no.na))[["greatest.pw2"]]

    # Compute delta
    delta <- length(x.no.na) - 2 ^ z

    # If delta < 20 , delete "delta" values randomly from the series
    p <- c()
    if (delta < 20) {
      if (delta > 0) {
        x.no.na <- x.no.na[-sample(1:length(x.no.na), delta)]
      }

      # Perform the test and return p.value for w.noise (Ho: WN)
      p <- genwwn.test(x.no.na)[["p.value"]]
    } else{
      #print("here 1")

      # Create 2 series of length 2^z:
      #   - from position 1 to 2^z
      #   - from position (n - 2^z) to n
      x.no.na.1 <- x.no.na[1:2 ^ z]
      n <- length(x.no.na)
      x.no.na.2 <- x.no.na[(n - 2 ^ z + 1):n]

      p1 <- genwwn.test(x.no.na.1)[["p.value"]]
      p2 <- genwwn.test(x.no.na.2)[["p.value"]]

      p <- p1
      if (p2 <= p1) {
        p <- p2
      }
    }
  } else{
    print("Error: WN test cannot applied, length of series < 16")
  }

  # Add the Box test
  p.box <- Box.test(x, box.test.lag)$p.value

  # Return the smallest of the 2 p-values
  p.ret <- p
  if (!is.na(p.box) &
      !is.na(p) &
      p.box <= p) {
    p.ret <- p.box
  }
  print(" ===> Selected p.box")

  if (is.na(p.ret)) {
    stop("WN test cannot applied, length of series < 16")
  }

  print(paste("WM test p-value=", p.ret, sep=""))
  return(p.ret)
}

###################
# smoothing spline
###################
apply.gam.hastie <- function(y, fun = "spline", k) {
  x <- 1:length(y)

  # Use spline or loess in gam
  fun.type <- "s"
  if (fun == "loess") {
    fun.type -> "lo"
  }
  #form <- paste("y ~ ", fun.type,"(x)", sep="") s(x, k = -1, bs = "cs"
  form <- paste("y ~ ", fun.type, "(x, ", k , ")", sep = "")
  m <- gam(as.formula(form), na.action = na.exclude)

  return(predict(m, newdata = data.frame(x = x)))
}



###########################
# Validate input dataset
###########################
validate.input <- function(df, group.name) {
  # Return value
  ret <- NA

  # Verify df is matrix or df
  is.matrix.or.df <- TRUE

  if (!is.matrix(df) & !is.data.frame(df)) {
    stop("\n- \"met.dat\" must be a matrix or a data frame")
    #quit(save="no")
  }

  # All Missing values ?
  mis <- unlist(lapply(df, function(x) {
    return(sum(is.na(x)))
  }))
  if (length(mis[mis == 0]) == 0) {
    stop("\n- missing values detected in ALL columns of \"met.dat\"")
  }

  flag.rownames                  <- FALSE
  flag.rownames.group.name       <- FALSE
  flag.rownames.group.number     <-
    FALSE # is.numeric("group numbers")?
  flag.rownames.group.number.seq <-
    FALSE # is.numeric("group numbers") AND in sequential order?
  flag.rownames.order.name       <- FALSE # is this == "order"
  flag.rownames.order.number     <-
    FALSE # is.numeric("order numbers")?
  flag.rownames.order.number.seq <-
    FALSE # is.numeric("order numbers") AND in sequential order?

  # Verify rownames structure
  z <- sapply(rownames(df), function(x) {
    return(strsplit(x, "_"))
  })

  # (1) verify that the length of each vector in the list z is
  # either 4 ("plate_<>_order_<>") or 6 ("plate_<>_order_<>_id_<>")
  z.l <- unlist(lapply(z, function(x) {
    return(length(x))
  }))
  if (!(sum(z.l) == 4 * length(z.l) | sum(z.l) == 6 * length(z.l))) {
    flag.rownames <- TRUE
  }

  #### If the row names are wrong STOP here and return an error
  # stopifnot(!flag.rownames)
  if (flag.rownames) {
    stop(
      "\n - wrong row names format must be: \"<group name>_<group number>_order_<order number>\" or \"<group name>_<group number>_order_<order number>_id_<participant id>\""
    )
  }

  # (2) Validate fields
  group.names     <- unlist(lapply(z, function(x) {
    return(x[1])
  }))
  group.numbers   <- unlist(lapply(z, function(x) {
    return(x[2])
  }))
  order.names     <- unlist(lapply(z, function(x) {
    return(x[3])
  }))
  order.numbers   <- unlist(lapply(z, function(x) {
    return(x[4])
  }))

  # Check group names are all the same
  if (sum(group.names == group.name) != length(group.names)) {
    flag.rownames.group.name <- T
  }

  # Check the "group numbers" are actually alphanumeric
  is.num <-
    sapply(group.numbers, function(x) {
      return(!grepl("^[A-Za-z]+$", x , perl = T))
    })
  if (sum(is.num) != length(group.numbers)) {
    flag.rownames.group.number <- T
  }
  #print(flag.rownames.group.number)

  # Check that the group numbers are in sequential order
  if (!flag.rownames.group.number) {
    if (is.unsorted(as.numeric(group.numbers))) {
      flag.rownames.group.number.seq <- T
    }
  } else{
    flag.rownames.group.number.seq <- T
  }

  # Check "order" names
  if (sum(order.names == "order") != length(order.names)) {
    flag.rownames.order.name <- T
  }

  # Check "order numbers" are alphanumeric
  is.num <-
    sapply(order.numbers, function(x) {
      return(!grepl("^[A-Za-z]+$", x , perl = T))
    })
  if (sum(is.num) != length(order.numbers)) {
    flag.rownames.order.number <- T
  }
  #print(flag.rownames.order.number)
  if (!flag.rownames.order.number) {
    if (is.unsorted(as.numeric(order.numbers))) {
      flag.rownames.order.number.seq <- T
    }
  } else{
    flag.rownames.order.number.seq <- T
  }


  # Check that the order numbers are in sequential order
  # if(is.unsorted(as.numeric(order.numbers))){flag.rownames.order.number.seq <- T}

  messages <- c()
  if (flag.rownames.group.name) {
    messages <-
      c(
        messages,
        paste(
          "- wrong group_name field in one or more row names (expected:\"",
          group.name,
          "\")",
          sep = ""
        )
      )
  }
  if (flag.rownames.group.number) {
    messages <-
      c(messages,
        "- group_number field not alphanumeric in one or more row names")
  }
  if (flag.rownames.group.number.seq) {
    messages <- c(messages, "- group_number(s) not sorted in row names")
  }

  if (flag.rownames.order.name) {
    messages <-
      c(messages,
        "- wrong order_name field in one or more row names (expected: \"order\")")
  }
  if (flag.rownames.order.number) {
    messages <-
      c(messages,
        "- order number field not alphanumeric in one or more row names")
  }
  if (flag.rownames.order.number.seq) {
    messages <-
      c(messages, "- order_number(s) not sorted in row names")
  }

  # Should we issue an error
  issue.error <-
    any(
      flag.rownames.group.name,
      flag.rownames.group.number,
      flag.rownames.group.number.seq,
      flag.rownames.order.name,
      flag.rownames.order.number,
      flag.rownames.order.number.seq
    )

  if (issue.error) {
    # Create stop message
    stop.message <- paste("\n", paste(messages, collapse = "\n"))
    stop(stop.message)
  }

  # If we made it so far, then we can return the dataset with complete
  ret <- df
  if (length(names(mis[mis == 0])) != dim(df)[2]) {
    complete <- names(mis[mis == 0])
    ret <- df[complete]
    print(
      paste(
        "Warning: found",
        length(complete),
        "complete (i.e. with no NA) variables out of",
        dim(df)[2]
      )
    )
  }
  return(ret)
}






#################################################################
# Returns named boolean vector
# - TRUE : failed to reject homogeneity of variance
# - FALSE: reject homogeneity of variance
#################################################################
homogen.var <- function(met.dat)
{
  # Create a vector to store results of homogeneity variance test
  is.hom.var <- rep(FALSE, dim(met.dat)[2])
  names(is.hom.var) <- names(met.dat)
  #names(is.hom.var) <- vars

  #str(met.dat)
  #print(rownames(met.dat))

  for (i in names(met.dat)) {
    #for(i in vars){
    met.df <-  met.dat[i]

    if(sd(met.dat[,i], na.rm=T) == 0){
      is.hom.var[i] <- TRUE}
    else{
      met.df[["group"]] <- sapply(rownames(met.df),
                                  function(x) {
                                    strsplit(x, "_")[[1]][2]
                                  })

      # Shapiro test of normality
      x <- shapiro.test(met.df[, i])

      # If normality  use Levine test for homogeneity of variance
      # Else us Fligner-Killeen test
      p.val <- c()
      if (x$p.value < 0.05) {
        # no normality : Fligner-Killeen
        flig.tst <-
          fligner.test(as.formula(paste(i, " ~ group", sep = "")), data = met.df)
        p.val <- flig.tst$p.value
      } else{
        # Change this once "car" library is available
        lev.tst <- levene.test(met.df[, i],
                               met.df[, "group"],
                               location = "median",
                               correction.method = "zero.correction")
        p.val <- lev.tst$p.value
      }
      if (p.val > 0.05) {
        is.hom.var[i] <- TRUE
      }
    }
  }

  return(is.hom.var)
}


################################################
# Residualize
################################################
residualize <- function(met.df, met, group.var, plates) {
  for (i in 1:(length(plates) - 1)) {
    on.plate <- paste("in.", i, sep = "")
    met.df[[on.plate]] <- 0
    met.df[grep(paste(group.var, "_", plates[i], "_", sep = ""),
                rownames(met.df)),][[on.plate]] <- 1
  }

  # Formula
  form <-
    paste(met, " ~ ", paste(paste(
      "as.factor(in.", 1:(length(plates) - 1), ")", sep = ""
    ), collapse = " + "))

  # Run regression and residualize if no errors
  ret <- tryCatch({
    reg <- lm(as.formula(form), data = met.df, na.action = na.exclude)

    # Calling "resid" on the regression object created with the na.action=na.exclude
    # option takes care of the missing values. "resid(reg)" return NA if the
    # observation is also NA (remember that reg$resid will only return the residuals for the
    # non-missing values)
    list(res = resid(reg), status = 0)
  },
  error = function(e) {
    print(paste("Error processing", met))
    print(e)

    # Residualization failed, return a vector of "res.failed"
    return(list(res = rep("res.failed", dim(met.df)[1]), status = -1))
  }) # End of try-catch
  return(ret)
}

######################################################
# Determine number of independent dimensions by PCA
######################################################
do.pca <- function(res.norm.met.dt) {
  # PCA
  n.90.pc <- tryCatch({
    pc <- prcomp(res.norm.met.dt, scale = TRUE)
    cum.prop    <- summary(pc)$importance["Cumulative Proportion",]
    print(paste("PCA - number of dimensions:", min(which(cum.prop >= 0.9))))
    min(which(cum.prop >= 0.9))
  }, warning = function(w) {
    print(paste("WARNING @ PCA:", w))
    return(NA)
  }, error = function(e) {
    print(paste("ERROR @ PCA", e))
    return(NA)
  })

  ret.n.90.pc <- c()
  if (is.na(n.90.pc)) {
    ret.n.90.pc <- dim(res.norm.met.dt)[2]
  }
  else{
    ret.n.90.pc <- n.90.pc
  }

  return(ret.n.90.pc)
}

#######################################
# Returns original name of metabolite
######################################
get.orig.met.name <- function(met) {
  met.original <- c()
  if(str_starts(string = met, pattern = "resid\\.norm\\.")) {
    met.original <-str_remove(string = met, pattern = "^resid\\.norm\\.")
  }else if(str_starts(string = met, pattern = "resid\\.")) {
    met.original <- str_remove(string = met, pattern = "^resid\\.")
  }else if(str_starts(string = met, pattern = "norm\\.")) {
    met.original <- str_remove(string = met, pattern = "^norm\\.")
  }else if(str_starts(string = met, pattern = "detr\\.resid\\.norm\\.")) {
    met.original <- str_remove(string = met, pattern = "^detr\\.resid\\.norm\\.")
  }else if(str_starts(string = met, pattern = "detr\\.resid\\.")) {
    met.original <- str_remove(string = met, pattern = "^detr\\.resid\\.")
  }else if(str_starts(string = met, pattern = "detr\\.norm\\.")) {
    met.original <- str_remove(string = met, pattern = "^detr\\.norm\\.")
  }else if(str_starts(string = met, pattern = "detr\\.")) {
    met.original <- str_remove(string = met, pattern = "^detr\\.")
  }else{
    met.original <- met
  }

  return(met.original)
}


##########################################
# Validate formula for acf guided spline
##########################################
validate.formula <- function(form, DEBUG=FALSE){
  valid          <- TRUE
  fmla.out       <- NULL
  fmla.spl.type  <- NULL
  fmla.spl.term  <- NULL
  fmla.expl      <- NULL
  fmla.covs      <- NULL

  f.str <- as.character(form)
  outcome <- f.str[2]
  other <- unlist(strsplit(gsub(" ", "", f.str[3], fixed = TRUE),"\\+"))

  w.covs <- FALSE
  if(length(other) > 1){w.covs <- TRUE}

  # Verify there is only one spline term in the formula ("s" or "bs")
  bss <- length(grep("^bs\\(", other, value=T))
  ss  <- length(grep("^s\\(", other, value=T))

  if(DEBUG){
    print(paste("bss=",bss))
    print(paste("ss=",ss))
    print(paste("valid=", valid))
  }

  # Issue error if formula is not valid
  if(bss == 0 & ss == 0){
    print("Error: no spline terms found in formula")
    valid <- FALSE
  }else if((bss != 0 & bss > 1) | (ss != 0 & ss > 1)){
    print("Error: only 1 spline term is allowed in formula");
    valid <- FALSE
  }else if(bss + ss > 1){
    print("Error: only 1 spline term is allowed in formula");
    valid <- FALSE
  }else{
    # Formula is ok
    if(bss == 1){
      fmla.spl.type <- "b.spline"
    }else{
      fmla.spl.type <- "s.spline"
    }

    # Spline term
    fmla.spl.term <- grep("s\\(", other, value=T)

    # Expl. variable
    fmla.expl <- strsplit(strsplit(fmla.spl.term,"\\(")[[1]][2], ")")[[1]][1]

    # Covariates
    if(w.covs){
      fmla.covs <- other[!other %in% fmla.spl.term]
    }
  }

  if(DEBUG){
    print(paste("valid (post-check):", valid))
    print(paste("fmla.spl.type:", fmla.spl.type))
    print(paste("fmla.spl.term:", fmla.spl.term))
    print(paste("fmla.expl:", fmla.expl))
    print(paste("w.covs:", w.covs))
    print(paste("fmla.covs:", fmla.covs))
  }

  return(list(fmla.valid=valid, fmla.out=outcome, fmla.spl.type=fmla.spl.type,
              fmla.spl.term=fmla.spl.term, fmla.expl=fmla.expl, fmla.covs=fmla.covs))
}

###########################
# ACF guided spline
###########################
atf.mgcv.spline <- function(formula, data, df.range, spl.type="cr", boxlag=20, DEBUG=FALSE){

  # Print function call
  print(match.call())

  # Validate formula
  fmla <- validate.formula(formula, DEBUG)

  # Exit if formula is not valid
  if(!fmla$fmla.valid){
    return(0)
  }

  # Get the expl. variable
  expl.str <- fmla$fmla.expl

  # Sort data frame by expl. variable
  data <- data[order(data[[expl.str]]), ]

  # Extract explanatory variable
  expl <- data[[expl.str]]

  # Extract covariates
  covs <- ""
  str(fmla$fmla.covs)
  if(!is.null(fmla$fmla.covs)){
    covs <- paste(" + ", paste(fmla$fmla.covs, collapse=" + "), sep="")
    if(DEBUG){
      print(paste("covs:", covs))
    }
  }

  # Spline basis type
  if(spl.type == "cr"){
    bs <- "bs='cr'"
  }else if(spl.type == "tp"){
    bs <- "bs='tp'"
  }

  # Set the fx term (fx=TRUE for unpenalized spline)
  fx <- TRUE

  # Loop through the dfs
  box.pvals <- c(); best.fit <- NULL; best.k <- c();
  resp <- c(); spl.term <- c()

  first.k <- TRUE
  for(k in df.range){
    spl <- paste("s(", expl.str,", ",bs,", k=",k, ", fx=",fx,")" , covs, sep="")
    new.f <- paste(fmla$fmla.out, "~", spl, sep=" ")

    if(DEBUG){
      print(paste("Formula: ", new.f))
    }

    # Run the regression
    fit <- mgcv::gam(as.formula(new.f), data=data)

    # Residuals
    res <- residuals(fit)

    # Ljung-Box test of autocorrelations on the residuals
    pvalue <- Box.test(res, lag=min(boxlag,length(res)-1), type = "Ljung-Box")$p.value

    # Is this the best fit ?
    if(first.k){
      best.fit <- fit
      best.k   <- k
      first.k <- FALSE
    }else{
      # Decide if strictly ">" or ">="
      if(pvalue >= max(box.pvals)){
        best.fit <- fit
        best.k   <- k
        print(paste("Best k:", best.k))
        print(paste("Best p:", pvalue))
      }
    }
    box.pvals <- c(box.pvals,pvalue)
  }

  # Add names to vector of pvalues
  names(box.pvals) <- df.range

  if(DEBUG){
    print(paste("Best param:", best.k))
    plot(expl, data[[fmla$fmla.out]], col="grey", xlab=expl.str, ylab=fmla$fmla.out)
    lines(expl, best.fit$fitted.values, col="red")
  }

  if(max(box.pvals) < 0.05){warning("P-value relative to best df < 0.05")}
  return(list(best.par=best.k, best.pval=max(box.pvals), box.pvals=box.pvals, best.fit=best.fit))
}

######################
# Normalize variance #
######################
normalize.variance <- function(met.dat, vars, summary.transf, group.var, final=FALSE){

  #########################################################
  # Step 1: Test for homogeneity of variance among plates #
  #########################################################
  print("=> Testing homogeneity of variance and normalizing ...")
  print(vars[1:10])
  met.dat <- met.dat[vars]
  str(met.dat)
  # plot(met.dat[[vars]])
  # print(sum(is.na(met.dat[[vars]])))
  plates <- unique(unlist(lapply(strsplit(rownames(met.dat), "_"), function(x) {
    return(x[2])
  })))
  if (length(plates)<=1){
    is.hom.var<-rep(TRUE,dim(met.dat)[2])
    names(is.hom.var) <- names(met.dat)
  }else{
    is.hom.var <- homogen.var(met.dat)
  }
  # Update status for metsbs that pass hom. var.
  mets.pass.hom.var <- names(is.hom.var[is.hom.var])
  print(mets.pass.hom.var)
  for(i in mets.pass.hom.var){
    if(final){
      summary.transf[summary.transf$metabolite == get.orig.met.name(i), "norm.detr"] <- FALSE
    }else{
      #print("HHHHHHHHHHH")
      summary.transf[summary.transf$metabolite == get.orig.met.name(i), "normalized"] <- FALSE
    }
  }

  #print(table(is.hom.var), exclude=NULL)

  ############################################################
  # Step 2: Normalize by plate SD if homogeneity of variance
  # fails
  ############################################################
  # Get metabs that do not pass hom. var. test
  mets.no.hom.var <- names(is.hom.var[!is.hom.var])

  print("mets.no.hom.var")
  str(mets.no.hom.var)
  print("===========")
  # Loop through the metabolites and normalize by sd of groups (e.g. "plates")
  plates <-
    unique(unlist(lapply(strsplit(
      rownames(met.dat), "_"
    ), function(x) {
      return(x[2])
    })))

  print("PLATES:")
  str(plates)

  for (i in mets.no.hom.var) {
    met.plates.sd <- c()
    met.plates.mean<-c()
    print(plates)
    for (j in plates) {
      print(paste("plate:", j))
      grep.plates <-
        grep(paste(group.var, "_", j, "_", sep = ""), row.names(met.dat))
      sd.plate    <-  sd(met.dat[grep.plates, i], na.rm = T)
      mean.plate <- mean(met.dat[grep.plates, i], na.rm = T)
      print(paste("SD PLATE:", sd.plate))

      if(!is.na(sd.plate) & sd.plate == 0) {
        sd.plate <- NA
      }

      met.plates.sd <-
        c(met.plates.sd, rep(sd.plate, length(grep.plates)))
      met.plates.mean<-c(met.plates.mean, rep(mean.plate, length(grep.plates)))
    }

    # Normalize
    met.dat[[paste("norm.", i, sep = "")]] <- met.dat[[i]] / met.plates.sd
    str(met.dat)
    print(i)
    plot(met.dat[[i]], type="l")
    plot(met.dat[[paste("norm.", i, sep = "")]], type="l")

    # Update summary
    print(summary.transf)
    if(final){
      print("FINAL == TRUE")
      summary.transf[summary.transf$metabolite == get.orig.met.name(i), "norm.detr"] <-TRUE
      print("****************")
      print(summary.transf)
      print("****************")
    }else{
      print("FINAL == FALSE")
      summary.transf[summary.transf$metabolite == get.orig.met.name(i), "normalized"] <-TRUE
      print("****************")
      print(summary.transf)
      print("****************")
    }
  }

  return(list(met.dat = met.dat, summary.transf=summary.transf))
}


########
# WiNN
########
#' @title WiNN in development
#'
#' @description A function for metabolite correction
#' @param input.dat input dataset as data frame or matrix
#' @param input.dat Study ID.
#' @param group.var Grouping variable. Defaults to "plate".
#' @param max.knots Maximum number of knots as % .
#' @param debug Debug mode. Defaults to TRUE.
#' @param runall Apply correction regrdless of white noise test.Defaults to FALSE.
#' @param save.corrected.data Save corrected dataset. Defaults to TRUE
#' @param save.pdfs Save pdfs of corrected/uncorrected data. Defaults to TRUE
#' @return A list (corrected.strict, corrected.lenient,corrected.resid, correction.summary)
#' @keywords WiNN
#' @import hwwntest
#' @import gam
#' @import mgcv
#' @import stringr
#' @import lawstat
#' @examples
#' winn()
winn <-
  function(input.dat,
           group.var = "plate",
           max.knots = 10,
           debug = T,
           runall = F) {

    #### FOR DEVELOPMENT ONLY ##############
    if (FALSE) {
      dir.pp <-
        "/an/vital/vital200/QC/correction_gam_2020/JUNE_2020_replace_outliers/7_APPLY_TO_OTHER_DSETS/OTHER_DSETS/CVD_CACO"
      input.dat      <-
        get(load(paste(
          dir.pp, "/cvdcaco_met_uncorrected.RData", sep = ""
        )))
      pp.dt               <-
        get(load(paste(dir.pp, "/pp_no_outl.RData", sep = "")))
      stdy.id <- "CVD_CACO_12152020"
      group.var = "plate"
      max.knots = 10
      debug = T
      runall = F
      save.corrected.data = T
      save.uncorrected.data = T
      save.pdfs = T
    }
    ###########################

    print("Running WiNN with following parameters:")
    print(paste("met.dat =", deparse(substitute(input.dat))))
    print(paste("max.knots =", max.knots))
    print(paste("debug =", debug))
    print(paste("runall =", runall))
    print("########################################")

    ################################
    # Validate input
    ################################
    print("=> Validating input data...")
    met.dat <-
      tryCatch(
        validate.input(df = input.dat, group.name = group.var),
        error = function(e) {
          #message("An error occurred:\n", e)
          message("", e)
          quit(save = "default")
        },
        warning = function(w) {
          message("A warning occured:\n", w)
        }
      )

    #################################
    # Create study directory
    #################################

    ########################################################
    # Create a summary data frame to keep track of
    # the transformations/corrections performed on the data
    ########################################################
    n.mets <- dim(met.dat)[2]
    summary.transf <- data.frame(
      metabolite = names(met.dat),
      normalized = rep(NA, n.mets),
      residualized = rep(NA, n.mets),
      pass.wn.1 = rep(NA, n.mets),
      detrended = rep(NA, n.mets),
      norm.detr = rep(NA, n.mets),
      pass.wn.2 = rep(NA, n.mets))

    ###########################################
    # Step 1: Test for homogeneity of variance
    # and normalize if it fails
    ###########################################
    norm.v <- normalize.variance(met.dat=met.dat, vars=names(met.dat),
                                 summary.transf=summary.transf, group.var=group.var)

    met.dat        <- norm.v$met.dat
    summary.transf <- norm.v$summary.transf

    print(summary.transf)

    #############################################
    # Step 2:  Residualize by plate number
    #############################################
    # First loop through the original metabolite names and determine
    # which one should go through the ANOVA test and be potentially
    # residualized. If the metabolite has been normalized it is
    # the norm.<original metab name>" that should be tested,
    # otherwise "<original metab name>":
    #
    # sd.normalized
    #       F             "<original metab name>"
    #       T             "norm.<original metab name>"

    # Get the plates
    plates <- unique(unlist(lapply(strsplit(rownames(met.dat), "_"), function(x) {
      return(x[2])
    })))

    print("=> Anova and residualization ...")
    metabs.to.anova.test <- c()
    for (i in as.character(summary.transf$metabolite)) {
      x <- summary.transf[as.character(summary.transf$metabolite) == i, ]
      if (x$normalized == FALSE) {
        metab.to.test <- i
      }
      if (x$normalized == TRUE) {
        metab.to.test <- paste("norm.", i, sep = "")
      }
      metabs.to.anova.test <- c(metabs.to.anova.test, metab.to.test)
    }

    # First let's determine the presence of a plate effect
    n.anova <- 0
    n.residualized <- 0

    for (met in metabs.to.anova.test) {
      # Extract the original name if this metabolite was normalized
      met.original <- met
      if (str_starts(string = met, pattern = "norm.")) {
        met.original <- str_remove(string = met, pattern = "^norm\\.")
      }
      if (length(plates)<=1){
        a.pval<-1
      }else{
        # Create the group variable for the ANOVA
        met.df            <- met.dat[met]
        met.df[["group"]] <-
          as.factor(paste(group.var, "_", unlist(lapply(strsplit(rownames(met.df), "_"),
                                                        function(x) {return(x[2])})), sep = ""))

        # ANOVA test
        a <- aov(as.formula(paste(met, "~ group", sep = "")), data = met.df)
        a.pval <- summary(a)[[1]][, "Pr(>F)"][1]
      }
      # Residualize if ANOVA test is significant
      if (a.pval < 0.05) {
        n.anova <- n.anova + 1
        residualized <- residualize(met.df, met, group.var, plates)

        # Was the residualization successfull, i.e. no errors?
        # If so, replace the metabolite value with  the residuals and
        # add the "resid" prefix to the name
        if (residualized$status == 0) {
          n.residualized <- n.residualized + 1
          met.dat[[paste("resid.", met, sep = "")]] <-residualized$res

          # Update summary
          summary.transf[summary.transf$metabolite == met.original, "residualized"] <-
            TRUE
        }else{
          summary.transf[summary.transf$metabolite == met.original, "residualized"] <-
            FALSE
        }
      }else{
        summary.transf[summary.transf$metabolite == met.original, "residualized"] <-
          FALSE
      }
    } ## Closes "for (met in metabs.to.anova.test)"

    print(paste(n.residualized, " metabolites have been residualized", sep=""))

    ##########################################################################
    # Step 4: WN test (first)
    ##########################################################################
    #
    # First loop through the original metabolite names and determine
    # which one should go through the WN test and be potentially
    # corrected. This is based on the transformation the metab
    # went through so far.
    #
    # sd.normilized residualized   metab
    #       F           F         "<original metab name>"
    #       F           T         "resid.<original metab name>"
    #       T           F         "norm.<original metab name>"
    #       T           T         "resid.norm.<original metab name>"
    #########################################################################
    metabs.to.wn.test <- c()
    for (i in as.character(summary.transf$metabolite)) {
      x <- summary.transf[as.character(summary.transf$metabolite) == i, ]
      if (x$normalized == F &
          x$residualized == F) {
        metab.to.test <- i
      }
      if (x$normalized == F &
          x$residualized == T) {
        metab.to.test <- paste("resid.", i, sep = "")
      }
      if (x$normalized == T &
          x$residualized == F) {
        metab.to.test <- paste("norm.", i, sep = "")
      }
      if (x$normalized == T &
          x$residualized == T) {
        metab.to.test <- paste("resid.norm.", i, sep = "")
      }
      metabs.to.wn.test <- c(metabs.to.wn.test, metab.to.test)
    } ## Closes "for(i in ...)"

    # Test for WN
    print("=> Testing for WN...")
    p.wn.tst <- sapply(metabs.to.wn.test, function(x) {
      ret <- NA; p <- NA
      p <- tryCatch({
        test.wn(met.dat[, x])
      }, warning = function(w) {
        print(paste("WARNING @ test.wn > ", x, ":", w))
      }, error = function(e) {
        print(paste("ERROR @ test.wn > ", x, ":", e))
        return(NA)
      })

      if (!is.na(p)) {
        ret <- p
      }
      return(ret)
    })

    names(p.wn.tst)  <- metabs.to.wn.test

    # Determine number of independent dimensions by PCA
    is.wn.before.corr <- c()

    # OD!!! ORIGINAL WiNN CODE
    # # The number of independent dimensions
    # print("=> Running first PCA ...")
    # n.90.pc.before.corr <-
    #   do.pca(res.norm.met.dt = met.dat[metabs.to.wn.test])
    # if (debug) {
    #   print(paste(
    #     "DBG 1 - number independent components:",
    #     n.90.pc.before.corr
    #   ))
    # }

    # OD!!! global Ljung box test is not used anymore
    # Other changes are on lines 2120, 2286 and 2296
    # The number of independent dimensions
    n.90.pc.before.corr <- 50
    if (length(names(metabs.to.wn.test))>50) {
      print("=> Running first PCA ...")
      n.90.pc.before.corr <-
        min(50, do.pca(res.norm.met.dt = met.dat[metabs.to.wn.test]))
      if (debug) {
        print(paste(
          "DBG 1 - number independent components:",
          n.90.pc.before.corr
        ))
      }
    }
    # Set the p-value threshold to Bonferroni:
    p.thr.before.corr <- 0.05 / n.90.pc.before.corr
    p.thr.before.corr <- 1 # OD!!! force to detrend all
    if (debug) {
      print(paste("DBG 2 - p-value threshold:",  p.thr.before.corr))
    }

    # Determine which metab is WN (Ho: WN)
    is.wn.before.corr <-
      ifelse(p.wn.tst <= p.thr.before.corr, FALSE, TRUE)
    names(is.wn.before.corr) <- metabs.to.wn.test

    # Select the metabolites that cannot be tested for WN because
    # errors in the test for white noise function (eg length < 16)
    mets.cannot.test <-
      names(is.wn.before.corr[is.na(is.wn.before.corr)])
    if (length(mets.cannot.test) > 0) {
      print("Metabolites that could not be tested for WN:", log.file = cannot.correct)
      for (i in mets.cannot.test) {
        print(i, log.file = cannot.correct)
      }
    }

    # Select the metabolites that do *not* pass the WN test and must be detrended
    mets.to.correct <-
      names(is.wn.before.corr[!is.na(is.wn.before.corr) &
                                is.wn.before.corr == FALSE])

    print(paste(
      "Number of metabolites not passing WN to be detrended:",
      length(mets.to.correct)
    ))

    # Perform a second PCA restricted to the metabs to be detrended
    print("=> Running second PCA ...")
    n.90.pc.after.corr <-
      do.pca(res.norm.met.dt = met.dat[mets.to.correct])
    if (debug) {
      print(paste("DBG 3 - number independent components among mets to correct:",
                  n.90.pc.after.corr
      )
      )
    }

    p.thr.after.corr <- 0.05 / n.90.pc.after.corr
    if (debug) {
      print(paste(
        "DBG 4 - p-value threshold for mets to correct:",
        p.thr.after.corr
      ))
    }

    # Update the transf. summary with the results of the WN test
    original.names <- c()
    for (i in mets.to.correct){
      met.original <- c()
      if (str_starts(string = i, pattern = "resid\\.norm\\.")) {
        met.original <- str_remove(string = i, pattern = "^resid\\.norm\\.")
      } else if (str_starts(string = i, pattern = "resid\\.")) {
        met.original <- str_remove(string = i, pattern = "^resid\\.")
      } else if (str_starts(string = i, pattern = "norm\\.")) {
        met.original <- str_remove(string = i, pattern = "^norm\\.")
      } else{
        met.original <- i
      }

      original.names <- c(original.names, met.original)

      summary.transf[summary.transf$metabolite %in% met.original, "pass.wn.1"] <-
        FALSE
    } ## Closes "for (i in mets.to.correct)"

    # Get  the metabs that pass WN test 1 and update summary transfer table
    print(original.names)
    pass.wn.1 <- as.character(summary.transf$metabolite[!summary.transf$metabolite %in% original.names])
    summary.transf[summary.transf$metabolite %in% pass.wn.1, "pass.wn.1"] <- TRUE
    summary.transf[summary.transf$metabolite %in% pass.wn.1, "detrended"] <- FALSE

    ####################################
    # Step 5: Detrend
    ####################################
    plates <- as.numeric(plates)

    # If this is true, just run them all regardless of WN test
    if (runall) {
      mets.to.correct <- metabs.to.wn.test
      p.thr.after.corr <- 0.05 / (length(metabs.to.wn.test))
      n.90.pc.after.corr <- length(metabs.to.wn.test)

      # Also, we are going to correct them all, so let's update
      summary.transf["detrended"] <- rep(TRUE, n.mets)
    }

    # Get plate coordinates. It will be used later for plotting
    plate.coord <- c()
    for (i in 1:length(plates)) {
      tmp.df <-
        met.dat[grep(paste(group.var, "_", i, "_order_", sep = ""),
                     rownames(met.dat),
                     value = T), ]
      plate.coord <-
        c(plate.coord, tail(sapply(rownames(tmp.df), function(x) {
          plate.order <-
            strsplit(x, "_id_")[[1]][1]
          return(strsplit(plate.order, "_order_")[[1]][2])
        }), n = 1))
    }

    # Get also the sequence order for later
    seq.order <- as.numeric(sapply(rownames(met.dat), function(x) {
      plate.order <-
        strsplit(x, "_id_")[[1]][1]
      return(strsplit(plate.order, "_order_")[[1]][2])
    }))

    ################################
    # main detrending loop
    ################################
    if (length(mets.to.correct) > 0) {
      print("=> Detrending ")

      # Create a data frame to store the detrending best # of knots and
      # p-values
      log.detrend <- data.frame(matrix(NA, nrow=length(mets.to.correct)*length(plates), ncol=2))
      names(log.detrend) <- c("met", "plate.bestk.bestp")

      # Initialize counters
      count <- 0
      count.met.plate <- 0

      print("The following metabolites will be detrended:")
      print(mets.to.correct)

      ###############################
      # Loop through the metabolites
      ###############################
      for (met in mets.to.correct) {

        # Get the original name
        met.original <- get.orig.met.name(met)

        # Extract the metabolite
        met.dt <- met.dat[met]

        print("##########################")
        print(paste("Detrending...", met))

        # Update count
        count <- count + 1
        if (count %% 20 == 0) {
          print(count)
        }

        pvals           <- c()
        p.best          <- 0

        # The vector that stores the "best" spline fits for each plate
        pred.trend <- c()

        # Loop through the plates
        for (i in plates) {
          count.met.plate <- count.met.plate + 1
          print(paste("Plate ", i))

          # Subset to this plate
          plate.dt.met <- met.dt[grep(paste(group.var, "_", i, "_", sep = ""),
                                      rownames(met.dt)), ]

          # OD!! add acf test to see if plate needs to be detrended
          plate.wn <- Box.test(plate.dt.met, lag=round(length(plate.dt.met)/2, digits=0), type = "Ljung-Box")$p.value
          print("!!!!###!!!! plate.wn=")
          print(plate.wn)


          # Run the knots optimization function
          mydf <- data.frame(x=1:length(plate.dt.met), y=plate.dt.met)

          # If the plate has only <=10 data points do *not* apply the spline regression
          # OD!!: if box ljung test for the plate is not signif do not run detrending either
          #if(dim(mydf)[1] <= 10 | sum(is.na(mydf$y)) > 0 ){
          #if(dim(mydf)[1] <= 10 | sum(is.na(mydf$y)) > 0 | plate.wn > .01){
          if(dim(mydf)[1] <= 10 | sum(is.na(mydf$y)) > 0 | plate.wn > .01){
            pred.trend <- c(pred.trend, rep(0, dim(mydf)[1]))
            log.detrend[count.met.plate, 1] <- met
            log.detrend[count.met.plate, 2] <- paste(i,":NA:NA")
          }else{
            n.stop <- (dim(mydf)[1]*max.knots) %/% 100 # <= the max number of knots
            print(paste("Dim plate: ", dim(mydf)[1], "; Knots: 3 - ", n.stop, sep=""))
            spl.fit <- atf.mgcv.spline(y ~ bs(x), data=mydf, 3:n.stop, DEBUG=TRUE)

            # append to pred.trend
            pred.trend <- c(pred.trend, spl.fit$best.fit$fitted.values)

            # Update the detrend log file
            log.detrend[count.met.plate, 1] <- met
            log.detrend[count.met.plate, 2] <- paste(i,spl.fit$best.par, spl.fit$best.pval, sep=":")
          }
        } # closes loop through plates

        # "Correct" the metabolite by subtracting the "best" spline fit from the metab
        z <- met.dt[[met]] - pred.trend

        # Define a new "detrended  metabolite" variable
        corr.met <- paste("detr.", met, sep = "")
        met.dt[[corr.met]] <- z

        # Add the detrended  metabolite to the original dset
        met.dat[[corr.met]] <- met.dt[[corr.met]]

        # Finally update status
        summary.transf[summary.transf$metabolite %in% met.original, "detrended"] <-
          TRUE

        ###################################################
        # Per Olga: test again for homog. variance and if
        # does not pass, normalize again
        ###################################################
        print(summary.transf)
        final.n.v <- normalize.variance(met.dat = met.dat,
                                        vars = corr.met,
                                        group.var=group.var,
                                        summary.transf = summary.transf,
                                        final = TRUE)
        norm2.corr.met <- NA
        if(dim(final.n.v$met.dat)[2] == 2){
          print("DBG 1")
          norm2.corr.met <- paste("norm", corr.met, sep=".")
          met.dt[[norm2.corr.met]] <- final.n.v$met.dat[[norm2.corr.met]]
        }

        summary.transf <- final.n.v$summary.transf

        str(summary.transf)
        print(summary.transf)

        # Get the white noise test p-value
        p.wn <- Box.test(z, lag=20, type = "Ljung-Box")$p.value

        # Now compare the "p.wn" value for this metabolite with the
        # threshold p.value determined in step 4
        result <- "Keep Ho: WN"
        if (p.wn <= p.thr.after.corr) {
          result <- "Reject Ho: not WN"
        }

        print(paste(corr.met, " :  p.wn=", p.wn, " - ", result, sep = ""))

        # Finally update status
        passwn2 <- p.wn > p.thr.after.corr
        summary.transf[summary.transf$metabolite %in% met.original, "pass.wn.2"] <- passwn2
        print(summary.transf)

        #################
        # Generate plots
        #################


      } # Closes "for(met in ...)"

      # Dump the detrend log to file


    } # Closes "if (length(mets.to.correct) > 0)"

    print(summary.transf)

    ############################################
    # Save transformed and corrected dataset
    ############################################
    # Loop through the metabolites and select the corresponding untransformed/transformed/corrected
    # columns according to the following criteria
    #
    #    corrected  sd.normalized   residualized    col name
    #       F           F               F           <original metab name>
    #       F           F               T           resid.<original metab name>
    #       F           T               F           norm.<original metab name>
    #       F           T               T           resid.norm.<original metab name>
    #       T           F               F           corr.<original metab name>
    #       T           F               T           corr.resid.<original metab name>
    #       T           T               F           corr.norm.<original metab name>
    #       T           T               T           corr.resid.norm.<original metab name>
    #

    col.names.sel <- c()
    corrected.df <- met.dat
    print(summary.transf$metabolite)

    for (i in summary.transf$metabolite) {
      x <- summary.transf[as.character(summary.transf$metabolite) == i, ]
      #x <- correction.summary[as.character(summary.transf$metabolite) == i, ]
      col.name <- c()
      if (x$detrended == F &
          x$normalized == F &
          x$residualized == F) {
        col.name <- i
      }
      if (x$detrended == F &
          x$normalized == F &
          x$residualized == T) {
        col.name <- paste("resid.", i, sep = "")
      }
      if (x$detrended == F &
          x$normalized == T &
          x$residualized == F) {
        col.name <- paste("norm.", i, sep = "")
      }
      if (x$detrended == F &
          x$normalized == T &
          x$residualized == T) {
        col.name <- paste("resid.norm.", i, sep = "")
      }

      if (x$detrended == T &
          x$normalized == F &
          x$residualized == F) {
        col.name <- paste("detr.", i, sep = "")
      }

      if (x$detrended == T &
          x$normalized == F &
          x$residualized == T) {
        col.name <- paste("detr.resid.", i, sep = "")
      }
      if (x$detrended == T &
          x$normalized == T &
          x$residualized == F) {
        col.name <- paste("detr.norm.", i, sep = "")
      }
      if (x$detrended == T &
          x$normalized == T &
          x$residualized == T) {
        col.name <- paste("detr.resid.norm.", i, sep = "")
      }

      corrected.df[i] <- met.dat[col.name]
      str(corrected.df)
      print("CHECKPOINT")

    } # Closes for (i in summary.transf$metabolite)

    corrected.df <- corrected.df[summary.transf$metabolite]
    print(names(corrected.df))


    # The dataset to save and return
    uncorrected.and.corrected <- met.dat



    ############################################
    # Save summary
    ############################################


    ##################################################
    # Create also pdfs of the metabolites that do not
    # go through detrending but are either:
    # - uncorrected
    # - sd-plate normalized
    # - residualized
    # - sd-normalized and residualized
    # dir.create(pdfs.d)
    # dir.create(untr.d)
    # dir.create(sd.norm.only.d)
    # dir.create(resid.only.d)
    # dir.create(sd.norm.resid.d)
    ##################################################


    # Returned object
    return(list(summary.transf = summary.transf, transf.corrected = corrected.df))
  }

