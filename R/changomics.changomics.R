
#####
#' @title changomics
#' @description A function for metabolite correction. This function first estimate the change points
#'  (jumps) and then correct each different segment by detrending and normalizing when needed.
#' @param data input dataset of metabolites as data frame.
#' @param max.knots Maximum number of knots as % .
#' @param debug Debug mode. Defaults to TRUE
#' @param runall Apply correction regrdless of white noise test.Defaults to FALSE.
#' @return the corrected dataframe of metabolite(s).
#' @import sarima
#' @import splines
#' @import dplyr
#' @importFrom stats na.omit
#' @export
#' @examples
#' met1<-rnorm(200)
#' changomics(as.data.frame(met1))

changomics<-function(data,
                     max.knots = 10,
                     debug = F,
                     runall = F){
  #check if data is a data frame. If not STOP algo.
  if (!is.data.frame(data)){
    stop("data must be a data frame")
  }
  NA_present<-FALSE
  #Check if there are NA's in the data frame.
  for (k.na in 1:dim(data)[2]) {
    if (sum(is.na(data[,k.na]))){
      NA_present<-TRUE
    }
  }
  #If there are Na's in the data frame, we omit the rows where there is at least one NA
  if (NA_present){
    warning("There is(are) NA value(s) in the data frame, na.omit() removed the rows containing NA(s).
            recommended-> remove Na's yourself to not remove a full row")
    data<-na.omit(data)
  }
  #If the data frame is empty after this step we stop the algo.
  if (nrow(data)==0){
    stop("All rows have been removed because there were Na's in each row: remove yourself Na's and run algo
         for one metabolites at a time")
  }
  corrected<-data # We create a second data frame where we will insert the corrected metabolites.
  for (k in 1:dim(data)[2]){
  knots<-seq(1,length(data[[k]]),length.out=floor(length(data[[k]])/60)+2) # We generate knots equally spaced. They will be used to estimate the segment moving mean
  knots<-knots[-c(1,length(knots))] # We take out first and last to get only the inner knots.
  knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots) # This step ensures that the knot will never be on the boundary (boundary will always be round number)
  changepoints<-fkPELT(data[[k]],knots = knots) # We estimate the change points using our developed algorithm
  colname<-colnames(data[k])
  if (is.null(colname)){colnames(data[k])<-paste("met_",k,sep="")}#If the column has no name we call it "met_met.number"
  df<-changetometa(data[k],changepoints) #We arrange the metabolite data frame to a data frame with
                                        #each row as "plate_number.plate_order_order.number" (this is to use winn algo properly)
  winnresult<-winn(df,group.var = "plate",
                   max.knots = max.knots,
                   debug = debug,
                   runall = runall) #With winn algo we correct the metabolites knowing their "plates/change points"
  corrected[k]<-winnresult$transf.corrected # We add the corrected result to the corrected data frame
  }
  return(corrected)

}
#########
# Change data to be used by winn algo
# This algo transform a data frame such that each row name is of the form "plate_plate.number_order_orde.number"
# data is data frame here
# Changepoints are the changepoints estimated in the sequence. They will delimit the different plates.
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
  ret<-as.data.frame(data)
  colnames(ret)<-colnam
  row.names(ret)<-names
  return(ret)
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
  #if size is too small, we do not fit a spline we just use the mean as an estimate.
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
  newknots<-knots[knots<index2& knots>index1] #we take only the knots that are inside both start index and end index
  spline<-lm(data~ns(cov,knots = newknots,intercept=TRUE)) # we estimate the moving mean (drift) with a spline
  mu<-spline$fitted.values
  sd<-sd(data-mu) #from this estimated moving mean we can estimate the sd
  # If sd is too small we return 0 cost
  if (sd<=1e-5){
    return(0)
  }
  neglog<-2*(sum((data-mu)^2/(2*sd^2))+size*log(sd*sqrt(2*pi))) #log likelihood of the segment supposing it is close to a
                                                                #normal distribution (usually the case with metabolomics data set)
  return(neglog)


}

########
#Estimating changepoints in time series (with drifts and jumps) with fixed knots predefined
#This algo is PELT algorithm but with a slight modification in the pruning and
#the cost of a segement as defined with the previous function "fksplinecost"
fkPELT<-function(data,knots){
  # Check if data is empty
  if (is.null(data)){stop("data cannot be NULL in fkPELT")}
  n<-length(data)
  f<-numeric(length = n+1)
  f[1]<- -3*log(300) # 3*log(n) is the penalty constant tha we use in the PELT algorithm
  cp<-rep(list(numeric()),n+1)
  R<-rep(list(numeric()),n)
  R[[1]]<-0


  for (t in 1:n){
    m<-length(R[[t]])
    neglog<-numeric(m)
    for (r in 1:m){
      neglog[r]<-fksplinecost(data[(R[[t]][r]+1):t],knots=knots,index1=R[[t]][r]+1,index2=t)
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

      tot<-f[R[[t]][r]+1]+log(300)+neglog[r]  #Here we have the pruning which is amplified by this "log(300)"
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


