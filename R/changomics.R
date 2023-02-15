#' @title changomics
#' @description A function for metabolite correction.
#' This function first estimate the change points
#'  (jumps) and then correct each different segment
#'  by detrending and normalizing when needed.
#' @param data input dataset of metabolites as data frame.
#' @param graph TRUE or FALSE. If TRUE will display the graph of uncorrected
#' data with estimated change points and correction
#' @return the corrected dataframe of metabolite(s).
#' @import splines
#' @import dplyr
#' @import lawstat
#' @importFrom graphics abline legend lines
#' @importFrom stats na.omit Box.test
#' @export
#' @examples
#' met1<-rnorm(200)
#' changomics(as.data.frame(met1))

changomics<-function(data,graph=FALSE){
  if (!is.data.frame(data)){
    stop("data must be a data frame") # We ensure the input is of type
                                      # data.frame
  }
  NA_present<-FALSE
  for (k.na in 1:dim(data)[2]) {
    if (sum(is.na(data[,k.na]))){
      NA_present<-TRUE
    }
  }
  if (NA_present){
    data<-na.omit(data) #remove NAs
    if (nrow(data)>0){
    warning("There is(are) NA value(s) in the data frame, na.omit()
    removed the rows containing NA(s).
            recommended-> remove Na's yourself to not remove a full row")
      }

  }
  if (nrow(data)==0){
    stop("All rows have been removed because there were Na's in each row:
    remove yourself Na's and run algo
         for one metabolites at a time") #If all values have been removed
    }
  corrected<-data
  for (k in 1:dim(data)[2]){
    knots<-seq(1,length(data[[k]]),length.out=floor(length(data[[k]])/60)+2)
    #Select number of knots ~1 knot per 60 observations
    knots<-knots[-c(1,length(knots))]
    knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
    changepoints<-fkPELT(data[[k]],knots = knots) # find the change points
    colname<-colnames(data[k])
    if (is.null(colname)){colnames(data[k])<-paste("met_",k,sep="")}
    #############
    #Detrending
    changes<-length(changepoints)
    datanew<-data[[k]]
    tau<-c(0,changepoints,length(data[[k]]))
    pred<-NULL #create a variable to store correction
    for (c in 1:(changes+1)){
      subset<-datanew[(tau[c]+1):tau[c+1]]
      size<-length(subset)
      if (size>=5){
        pvalue<-Box.test(subset,lag=floor(length(subset)/2),
                         type="Ljung-Box")$p.value
        # Test if the segment has autocorrelation
        if (!is.na(pvalue) & pvalue<0.05){
        dfw<-dfwhitenoise(subset) #Find the best degree of freedom
        spline<-smooth.spline(subset,df=dfw)
        pred<-c(pred,spline$y)#correction that will be applied
      }else{
        pred<-c(pred,rep(mean(subset),size)) #store the correction
      }

      }else{
        pred<-c(pred,rep(mean(subset),size))
      }
    }
    if (k==1){
    if (graph==TRUE){
      ## here is the code to plot visualizatio of how the algorithm works
      ## We display change points together with uncorrected signal and
      ## the correction that will be applied
      concentration<-data[[k]]
      plot(concentration,col="darkgray",
           xlab="reading sequence",ylab="uncorrected signal",
           main="uncorrected signal + correction and change points",
           type="l")
      for (l in changepoints) {
        abline(v=l,col="red",lty=2)
      }
      lines(1:length(data[[k]]),pred,col="red")
      legend(
        "bottomleft",
        lty = c(1, 1),
        col = c("darkgrey", "red"),
        legend = c("Uncorrected", "Correction"),
        cex=0.5
      )
    }}
    datanew<-datanew-pred
    #######
    #normalizing
    datanew<-changetometa(as.data.frame(datanew),changepoints)
    ## Put in form that we can know which segment/plate the data is in
    colnames(datanew)<-"orig"
    if (changes>=1){
      ## if no change no normalization
      homo2<-homogen.var.2(datanew) # Test if the different segments
                                    # have homogeneity of variance
      if (!homo2){

        plates <-
          unique(unlist(lapply(strsplit(
            rownames(datanew), "_"
          ), function(x) {
            return(x[2])
          })))

        met.plates.sd <- c()

        for (j in plates) {

          grep.plates <-
            grep(paste("plate", "_", j, "_", sep = ""), row.names(datanew))
          sd.plate    <-  sd(datanew[grep.plates,], na.rm = TRUE)



          if(!is.na(sd.plate) & sd.plate == 0) {
            sd.plate <- NA
          }

          met.plates.sd <-
            c(met.plates.sd, rep(sd.plate, length(grep.plates)))
        }

        # Normalize
        datanew[["norm"]] <- datanew[["orig"]] / met.plates.sd
        if (sum(is.na(datanew[["norm"]]))>0){

          dataret<-datanew[["orig"]]

        }

        dataret<-datanew[["norm"]]
      }
    }

    dataret<-datanew[["orig"]]


  #######
  #return result
    corrected[[k]]<-dataret
  }

  return(corrected)

}
#########
##change data to be in the form with group/plate to then normalize
changetometa<-function(data,changepoints){
  names<-vector("character",length = length(data[[1]]))
  tau<-c(0,changepoints,length(data[[1]]))
  for (i in 1:(length(tau)-1)){
    for (p in (tau[i]+1):tau[i+1]){
      plate_number<-i
      order<-p

      names[p]<-paste("plate",plate_number,"order",order,sep = "_")
    }
  }
  colnam<-colnames(data)
  ret<-as.data.frame(data[[1]])
  colnames(ret)<-colnam
  row.names(ret)<-names
  return(ret)
}

#############
#'@noRd
#'@importFrom stats smooth.spline

## Find the best degree of freedom that maximize the p-value of
## autocorrelation test
dfwhitenoise<-function(data){
  t<-n_distinct(data)
  l<-min(c(t-1,floor(t*20/100)))
  value<-numeric(length = l-1)
  for (i in 2:l){
    spline<-smooth.spline(data,df=i)
    epsilon<-data-spline$y
    value[i-1]<-Box.test(epsilon,lag=15,type="Ljung-Box")$p.value
  }
  df<-which.max(value)+1
  return(df)
}
###########
#' @noRd
#' @importFrom stats sd shapiro.test fligner.test as.formula

## Test if the different segments in the data set have the same variance
homogen.var.2 <- function(met.dat)
{
  # Create a vector to store results of homogeneity variance test
  is.hom.var <- FALSE

    met.df <-  met.dat["orig"]

    if(sd(met.dat[,1], na.rm=TRUE) == 0){
      is.hom.var <- TRUE}
    else{
      met.df[["group"]] <- sapply(rownames(met.df),
                                  function(x) {
                                    strsplit(x, "_")[[1]][2]
                                  })

      # Shapiro test of normality
      x <- shapiro.test(met.df[, 1])

      # If normality  use Levine test for homogeneity of variance
      # Else us Fligner-Killeen test
      p.val <- c()
      if (x$p.value < 0.05) {
        # no normality : Fligner-Killeen
        flig.tst <-
          fligner.test(as.formula(paste("orig", " ~ group", sep = "")),
                       data = met.df)
        p.val <- flig.tst$p.value
      } else{
        # Change this once bartlett.test from "car" library is available
        lev.tst <- levene.test(met.df[, 1],
                               met.df[, "group"],
                               location = "median",
                               correction.method = "zero.correction")
        p.val <- lev.tst$p.value
      }
      if (p.val > 0.05) {
        is.hom.var<- TRUE
      }
    }


  return(is.hom.var)
}

########
#Cost of segment calculator with fixed knots
#' @noRd
#' @importFrom stats lm
fksplinecost<-function(data,knots,index1=1,index2=length(data)){
  size<-length(data)
  if (size==1){
    return(0)
  }
  #if size is too small, we do not fit a spline
  #we just use the mean as an estimate.
  if (size<5){
    sd<-sd(data)
    mu<-mean(data)
    neglog<-2*(sum((data-mu)^2/(2*sd^2))+size*log(sd*sqrt(2*pi)))
    # If sd is too small we return 0 cost
    if (sd<=1e-5){
      return(0)
    }else{
      return(neglog)
    }
  }
  cov<-index1:index2
  newknots<-knots[knots<index2& knots>index1]
  # We take only the knots that are inside both start index and end index
  spline<-lm(data~ns(cov,knots = newknots,intercept=TRUE))
  # We estimate the moving mean (drift) with a spline
  mu<-spline$fitted.values
  sd<-sd(data-mu) #from this estimated moving mean we can estimate the sd
  # If sd is too small we return 0 cost
  if (sd<=1e-5){
    return(0)
  }
  neglog<-2*(sum((data-mu)^2/(2*sd^2))+size*log(sd*sqrt(2*pi)))
  #log likelihood of the segment supposing it is close to a
  #normal distribution (usually the case with metabolomics data set)
  return(neglog)


}



########
#Estimating changepoints in time series (with drifts and jumps)
# with fixed knots predefined
#This algo is PELT algorithm but with a slight modification in the pruning and
#the cost of a segement as defined with the previous function "fksplinecost"
fkPELT<-function(data,knots){
  # Check if data is empty
  if (is.null(data)){stop("data cannot be NULL in fkPELT")}
  n<-length(data)
  f<-numeric(length = n+1)
  f[1]<- -3*log(300) # 3*log(300) is the penalty constant tha we use
                     #in the PELT algorithm. Inspired by BIC criterion but
                     #modified
  cp<-rep(list(numeric()),n+1)
  R<-rep(list(numeric()),n)
  R[[1]]<-0


  for (t in 1:n){
    m<-length(R[[t]])
    neglog<-numeric(m)
    for (r in 1:m){
      neglog[r]<-fksplinecost(data[(R[[t]][r]+1):t],
                              knots=knots,index1=R[[t]][r]+1,index2=t)
    }
    stat<-numeric(m)
    for (r in 1:m){
      # Calculate the global cost with each possible change points in R[[t]]
      tot<-f[R[[t]][r]+1]+3*log(300)+neglog[r]
      stat[r]<-tot

    }
    f[t+1]<-min(stat)
    t1<-R[[t]][which.min(stat)]
    cp[[t+1]]<-c(cp[[t1+1]],t1)
    for (r in 1:m){

      tot<-f[R[[t]][r]+1]+log(300)+neglog[r]
      #Here we have the pruning which is amplified by this "log(300)"
      #, this makes the method a lot faster but a little bit less accurate
      if (tot<=f[t+1]& t<n){

        R[[t+1]]<-c(R[[t+1]],R[[t]][r])
      }
    }
    if (t<n&t>39){
      R[[t+1]]<-c(R[[t+1]],t-19)}
  }

  cp[[n+1]]<-cp[[n+1]][-1]
  return(cp[[n+1]]) #Return the set of change points
}







