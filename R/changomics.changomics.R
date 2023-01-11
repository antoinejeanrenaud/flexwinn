########
# PACKAGES to use
library(sarima)
library(dplyr)
library(splines)


#####
#' @title changomics
#' @description A function for metabolite correction
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
#' met1<-rnorm(1000)
#' changomics(as.data.frame(met1))

changomics<-function(data,
                     max.knots = 10,
                     debug = T,
                     runall = F){
  if (!is.data.frame(data)){
    stop("data must be a data frame")
  }
  NA_present<-FALSE
  for (k.na in 1:dim(data)[2]) {
    if (sum(is.na(data[,k.na]))){
      NA_present<-TRUE
    }
  }
  if (NA_present){
    warning("There is(are) NA value(s) in the data frame, na.omit() removed the rows containing NA(s).
            recommended-> remove Na's yourself to not remove a full row")
    data<-na.omit(data)
  }
  if (nrow(data)==0){
    stop("All rows have been removed because there were Na's in each row: remove yourself Na's and run algo
         for one metabolites at a time")
  }
  corrected<-data
  for (k in 1:dim(data)[2]){
  knots<-seq(1,length(data[[k]]),length.out=floor(length(data[[k]])/60)+2)
  knots<-knots[-c(1,length(knots))]
  knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
  changepoints<-fkPELT(data[[k]],knots = knots)
  colname<-colnames(data[k])
  if (is.null(colname)){colnames(data[k])<-paste("met_",k,sep="")}
  df<-changetometa(data[k],changepoints)
  winnresult<-winn(df,group.var = "plate",
                   max.knots = max.knots,
                   debug = debug,
                   runall = runall)
  corrected[k]<-winnresult$transf.corrected
  }

  return(corrected)

}
#########
#change data to be used by winn
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


##########
#' @noRd
#' @importFrom stats runif rbeta
datagen<-function(change=5){
  data<-numeric(0)
  datat<-numeric(0)
  datan<-numeric(0)
  v<-numeric(0)
  for (i in 1:(change+1)){
    mean<-runif(1,min=-1,max=1)
    max<-runif(1,1,2)
    trend<-sample(c(1,-1),1)
    sdfac<-runif(1,min=0.6,max=1)
    n<-sample(50:120,1)
    m<-sample(c(1,2,3,4),1)
    sdbeta<-sqrt((2*5)/((2+5)^2*(2+5+1)))
    if (m==1){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rnorm(n)#rt(n,10)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+trend*seq(0,max,length.out=n))
      datan<-c(datan,new+mean+trend*seq(0,max,length.out=n))

    }
    if (m==2){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+trend*seq(0,sqrt(max),length.out=n)^2)
      datan<-c(datan,new+mean+trend*seq(0,sqrt(max),length.out=n)^2)
    }
    if (m==3){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+trend*sqrt(seq(0,max,length.out=n)))
      datan<-c(datan,new+mean+trend*sqrt(seq(0,max,length.out=n)))

    }
    if (m==4){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean)
      datan<-c(datan,new+mean)

    }
    if (m==5){
      new<-(rbeta(n,2,5)-2/7)/sdbeta#rt(n,10)#rnorm(n)
      datat<-c(datat,new)
      data<-c(data,sdfac*new+mean+max*sin(seq(0,1,length.out=n)))
      datan<-c(datan,new+mean+max*sin(seq(0,10,length.out=n)))

    }
    v<-c(v,n)

  }
  v<-cumsum(v)
  v<-v[-length(v)]
  return(list(data,datat,v,datan))
}

######
#data  generation with no trend
#' @noRd
#' @importFrom stats runif rbeta
datagen.no.trend<-function(change=5){
  data<-numeric(0)
  datat<-numeric(0)
  datan<-numeric(0)
  v<-numeric(0)
  for (i in 1:(change+1)){
    mean<-runif(1,min=-4,max=4)
    sdfac<-runif(1,min=0.6,max=1)
    n<-sample(50:115,1)

    new<-(rbeta(n,2,5)-2/7)/0.16#rt(n,10)
    datat<-c(datat,new)
    data<-c(data,sdfac*new+mean)
    datan<-c(datan,new+mean)



    v<-c(v,n)
  }



  v<-cumsum(v)
  v<-v[-length(v)]
  return(list(data,datat,v,datan))
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
  if (size<5){
    sd<-sd(data)
    mu<-mean(data)
    neglog<-2*(sum((data-mu)^2/(2*sd^2))+size*log(sd*sqrt(2*pi)))
    if (sd==0){
      return(0)
    }else{
      return(neglog)
    }
  }
  cov<-index1:index2
  newknots<-knots[knots<index2& knots>index1]
  spline<-lm(data~ns(cov,knots = newknots,intercept=TRUE))
  mu<-spline$fitted.values
  sd<-sd(data-mu)
  if (sd==0){
    return(0)
  }
  neglog<-2*(sum((data-mu)^2/(2*sd^2))+size*log(sd*sqrt(2*pi)))
  if (size<15){
    return(neglog) }else{
      return(neglog)
    }

}

########
#Estimating changepoints in time series with fixed knots
fkPELT<-function(data,knots){
  n<-length(data)
  f<-numeric(length = n+1)
  f[1]<- -3*log(300)
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
      tot<-f[R[[t]][r]+1]+3*log(300)+neglog[r]
      stat[r]<-tot

    }
    f[t+1]<-min(stat)
    t1<-R[[t]][which.min(stat)]
    cp[[t+1]]<-c(cp[[t1+1]],t1)
    for (r in 1:m){

      tot<-f[R[[t]][r]+1]+log(300)+neglog[r]
      if (tot<=f[t+1]& t<n){

        R[[t+1]]<-c(R[[t+1]],R[[t]][r])
      }
    }
    if (t<n&t>39){
      R[[t+1]]<-c(R[[t+1]],t-19)}
  }

  cp[[n+1]]<-cp[[n+1]][-1]
  return(cp[[n+1]])
}

########
# Normalize when you know all changepoints correctly already




########
# Test if find correctly changepoints
#' @noRd
#' @importFrom graphics abline par
#' @importFrom stats ts
test<-function(nchange=6){
  x<-1:200
  y<-1:150
  #data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(146,mean=4,sd=1),rnorm(500,mean = sin(0.2*seq(1,125,by=0.25)),sd=0.5),rnorm(300,mean=seq(0.01,3,by=0.01)))
  #data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(200,mean = 0.0003*x^2,sd=1),rnorm(300,mean=0.5,sd=1.1),rnorm(300,mean = -0.5,sd=1.5),rnorm(100,mean = 0.5,sd=0.67))
  #data<-c(rnorm(450,mean = 0,sd=1),rnorm(150,mean=0.5,sd=1.1),rnorm(300,mean = 0,sd=1.5),rnorm(100,mean = 0.5,sd=1))
  #data<-rnorm(1000)
  data<-datagen(nchange)
  par(mfrow=c(2,1))
  plot(ts(data[[2]]),main="data untransformed")
  plot(ts(data[[1]]),main="data transformed with estimated changepoints")
  knots<-seq(1,length(data[[1]]),length.out=floor(length(data[[1]])/60)+2)
  knots<-knots[-c(1,length(knots))]
  knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
  changepoints<-fkPELT(data[[1]],knots = knots)
  for (b in data[[3]]){

    abline(v=b,col="blue")

  }
  correct<-0
  for (a in changepoints){
    if (sum(a==data[[3]])>0){
      correct<-1+correct
      abline(v=a,col="red")}else{
        abline(v=a,col="green")
      }
  }

  print(changepoints)
  print(correct/length(data[[3]]))
  par(mfrow=c(1,1))
}

