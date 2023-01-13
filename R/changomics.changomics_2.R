
#####
#' @title changomics.2
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

changomics.2<-function(data,
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
    # df<-changetometa(data[k],changepoints)
    # winnresult<-winn(df,group.var = "plate",
    #                  max.knots = max.knots,
    #                  debug = debug,
    #                  runall = runall)
    ch.result<-normalize(data[[k]],changepoints)
    corrected[[k]]<-ch.result
  }

  return(corrected)

}
#########
#change data to be used by winn
changetometa.2<-function(data,changepoints){
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
normalize<-function(data,changepoints){
  changes<-length(changepoints)
  datanew<-data
  tau<-c(0,changepoints,length(data))
  for (c in 1:(changes+1)){
    subset<-datanew[(tau[c]+1):tau[c+1]]
    size<-length(subset)

    if (size>=5){
      dfw<-dfwhitenoise(subset)
      spline<-smooth.spline(subset,df=dfw)
      datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]-spline$y
    }else{
      datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]-mean(datanew[(tau[c]+1):tau[c+1]])
    }



  }
  datanew<-changetometa.2(as.data.frame(datanew),changepoints)
  colnames(datanew)<-"orig"
  if (changes>=1){

    homo2<-homogen.var(datanew)
    if (!homo2){

      plates <-
        unique(unlist(lapply(strsplit(
          rownames(datanew), "_"
        ), function(x) {
          return(x[2])
        })))

      met.plates.sd <- c()
      print(plates)
      for (j in plates) {
        print(paste("plate:", j))
        grep.plates <-
          grep(paste("plate", "_", j, "_", sep = ""), row.names(datanew))
        sd.plate    <-  sd(datanew[grep.plates,], na.rm = T)

        print(paste("SD PLATE:", sd.plate))

        if(!is.na(sd.plate) & sd.plate == 0) {
          sd.plate <- NA
        }

        met.plates.sd <-
          c(met.plates.sd, rep(sd.plate, length(grep.plates)))
      }

      # Normalize
      datanew[["norm."]] <- datanew[["orig"]] / met.plates.sd
      if (sum(is.na(datanew[["norm."]]))>0){

        #p.wn <- Box.test(datanew[["orig"]], lag=20, type = "Ljung-Box")$p.value
        #pass.wn<- p.wn>0.05
        return(datanew[["orig"]])

      }

      #p.wn <- Box.test(datanew[["norm."]], lag=20, type = "Ljung-Box")$p.value
      #pass.wn<- p.wn>0.05

      return(datanew[["norm."]])
    }
  }
  #p.wn <- Box.test(datanew[["orig"]], lag=20, type = "Ljung-Box")$p.value
  #pass.wn<- p.wn>0.05
  return(datanew[["orig"]])

}
#############
#'@noRd
#'@importFrom stats smooth.spline
dfwhitenoise<-function(data){
  t<-n_distinct(data)
  l<-min(c(t-1,40))
  value<-numeric(length = l-1)
  for (i in 2:l){
    spline<-smooth.spline(data,df=i)
    epsilon<-data-spline$y
    auto<-autocorrelations(epsilon)
    auto.iid<-whiteNoiseTest(auto,h0="iid",x=epsilon)
    value[i-1]<-auto.iid[["test"]][3]
  }
  df<-which.max(value)+1
  return(df)
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







