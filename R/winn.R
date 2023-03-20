#' @title winn
#' @description A function for metabolite correction.
#' This function corrects batch effects and drifts
#' with the plate information by detrending and normalizing when needed.
#' @param data input dataset of metabolites as data frame.
#' @param graph TRUE or FALSE. If TRUE will display the graph of uncorrected
#' data with plates and correction and the corrected signal on a
#' second plot.
#' @param end.plates vector containing the position of the changes in plates
#' in terms of run order. The vector must be ordered.
#' @return the corrected dataframe of metabolite(s).
#' @import splines
#' @import dplyr
#' @import lawstat
#' @import mgcv
#' @importFrom graphics abline legend lines
#' @importFrom stats na.omit Box.test
#' @export
#' @examples
#' met1<-rnorm(200)
#' winn(as.data.frame(met1),end.plates=100)

winn<-function(data,end.plates,graph=FALSE){
  if (!is.data.frame(data)){
    stop("data must be a data frame") # We ensure the input is of type
    # data.frame
  }
  if (is.null(end.plates)){
    stop("end.plates cannot be empty") # We ensure there is at least two
    # plates
  }
  NA_present<-FALSE
  for (k.na in 1:dim(data)[2]) {
    if (sum(is.na(data[,k.na]))>0){
      NA_present<-TRUE
    }
  }
  if (NA_present){
      stop("There is(are) NA value(s) in the data frame, please remove Na's")

  }
  corrected<-data
  for (k in 1:dim(data)[2]){

    changepoints<-end.plates
    colname<-colnames(data[k])
    if (is.null(colname)){colnames(data[k])<-paste("met_",k,sep="")}
    changes<-length(changepoints)
    datanew<-data[[k]]-mean(data[[k]])
    #####
    #First normalization
    #This normalization is for the next anova test
    data.norm<-normalize.var(datanew,changepoints)
    ######
    #Residualizing if needed
    data.resid<-residualize(data.norm,changepoints)
    #############
    #Detrending
    tau<-c(0,changepoints,length(data.resid))
    pred<-NULL #create a variable to store correction
    for (c in 1:(changes+1)){
      subset<-data.resid[(tau[c]+1):tau[c+1]]
      size<-length(subset)
      if (size>=5){
        pvalue<-Box.test(subset,lag=floor(length(subset)/2),
                         type="Ljung-Box")$p.value
        # Test if the segment has autocorrelation
        if (!is.na(pvalue) & pvalue<0.01){
          dfw<-dfwhitenoise(subset)#Find the best degree of freedom
          x<-1:length(subset)
          spline<-mgcv::gam(subset~s(x,bs="cr",k=dfw,fx=TRUE))
          pred<-c(pred,spline$fitted.values)#correction that will be applied
        }else{
          pred<-c(pred,rep(0,size)) #store the correction
        }

      }else{
        pred<-c(pred,rep(0,size))
      }
    }
    if (k==1){
      if (graph==TRUE){
        ## here is the code to plot visualizatio of how the algorithm works
        ## We display change points together with uncorrected signal and
        ## the correction that will be applied
        concentration<-data.resid
        plot(concentration,col="darkgray",
             xlab="reading sequence",ylab="residualized/normalized signal",
             main="residualized/normalized signal + correction and plates",
             type="l")
        for (l in changepoints) {
          abline(v=l,col="red",lty=2)
        }
        lines(1:length(data.resid),pred,col="red")
        legend(
          "bottomleft",
          lty = c(1, 1),
          col = c("darkgrey", "red"),
          legend = c("residualized/normalized", "Correction"),
          cex=0.5
        )
      }}
    datanew<-data.resid-pred
    #######
    #normalizing by segment again
    data.norm.2<-normalize.var(datanew,changepoints)
    ######
    #Residualizing if needed
    data.resid.2<-residualize(data.norm.2,changepoints)
    #####
    #plot corrected
    if (k==1){
      if (graph==TRUE){

        ##We display the corrected signal for the first metabolite
        concentration<-data.resid.2
        plot(concentration,col="darkgray",
             xlab="reading sequence",ylab="corrected signal",
             main="corrected signal",
             type="l")

        legend(
          "bottomleft",
          lty = c(1, 1),
          col = c("darkgrey"),
          legend = c("corrected"),
          cex=0.5
        )
      }
    }
    #######
    #return result
    corrected[[k]]<-data.resid.2

  }

  return(corrected)

}

#'@noRd
#'@importFrom stats aov na.exclude
residualize<-function(met.data,endpoints){
  resid<-FALSE
  met.data<-changetometa(as.data.frame(met.data),endpoints)
  colnames(met.data)<-"orig"
  plates <-
    unique(unlist(lapply(strsplit(
      rownames(met.data), "_"
    ), function(x) {
      return(x[2])
    })))
    if (length(plates)<=1){
      a.pval<-1
    }else{
      met.data[["group"]] <-
        as.factor(paste("plate", "_",
                        unlist(lapply(strsplit(rownames(met.data), "_"),
                          function(x) {return(x[2])})), sep = ""))

      # ANOVA test
      a <- aov(as.formula(paste("orig", "~ group", sep = "")), data = met.data)
      a.pval <- summary(a)[[1]][, "Pr(>F)"][1]
    }
      if (a.pval<0.05){
        resid<-TRUE
      }


    if (resid==TRUE){
      for (i in 1:(length(plates) - 1)) {
        on.plate <- paste("in.", i, sep = "")
        met.data[[on.plate]] <- 0
        met.data[grep(paste("plate", "_", plates[i], "_", sep = ""),
                    rownames(met.data)),][[on.plate]] <- 1
      }

      # Formula
      form <-
        paste("orig", " ~ ", paste(paste(
          "as.factor(in.", 1:(length(plates) - 1), ")", sep = ""
        ), collapse = " + "))

      # Run regression and residualize if no errors
      ret <- tryCatch({
        reg <- lm(as.formula(form), data = met.data, na.action = na.exclude)

# Calling "resid" on the regression object created with the na.action=na.exclude
# option takes care of the missing values. "resid(reg)" return NA if the
# observation is also NA (remember that reg$resid will only return the
# residuals for the non-missing values)
        list(res = resid(reg), status = 0)
      },
      error = function(e) {
        print(paste("Error processing", "orig"))
        print(e)

        # Residualization failed, return a vector of "res.failed"
        return(list(res = rep("res.failed", dim(met.data)[1]), status = -1))
      }) # End of try-catch






      if (ret$status == 0) {
        data.resid<-ret$res

      }else{
        data.resid<-met.data[["orig"]]
      }
    }else{
      data.resid<-met.data[["orig"]]
    }


  return(data.resid)
}
