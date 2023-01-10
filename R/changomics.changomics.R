########
# PACKAGES to use
#library(tictoc)
library(sarima)
library(dplyr)
library(splines)
library(sigmoid)
library(sarima)
library(dplyr)
library(ggplot2)
library(ggpubr)
#####
#first modification of algo
changomics<-function(data,
                     max.knots = 10,
                     debug = T,
                     runall = F){
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

###########
#find best df for whitenoise

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

##########
#data generation
datagen.2<-function(size=800){
  data<-numeric(0)
  datan<-numeric(0)
  ch<-runif(1,min=floor(size/2)-100,max=floor(size/2)+100)
  sdbeta<-sqrt((2*5)/((2+5)^2*(2+5+1)))
  mean<-runif(2,min=-2,max=2)
  print(mean)
  max<-runif(2,1,2)
  trend<-sample(c(1,-1),2,replace = TRUE)
  sdfac<-runif(2,min=0.6,max=1)
  new1<-(rbeta(ch,2,5)-2/7)/sdbeta
  new2<-(rbeta(size-ch+1,2,5)-2/7)/sdbeta
  data<-c(new1,new2)
  datan<-c(sdfac[1]*new1+mean[1],sdfac[2]*new2+mean[2])
  return(list(datan,data,ch))
}

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
  # col<-ncol(ns(cov,knots = newknots,intercept=TRUE))
  # if (col>4) print(col)
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
#Estimating changepoints in distribution with fixed knots

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
    newR<-R[[t]]
    for (r in 1:m){
      neglog[r]<-fksplinecost(data[(newR[r]+1):t],knots=knots,index1=newR[r]+1,index2=t)
    }
    newneglog<-neglog
    # if (m<=10){
    #   newR<-R[[t]][R[[t]]<(t-(m-1))]
    #   newneglog<-neglog[R[[t]]<(t-(m-1))]
    #
    # }
    newm<-length(newR)
    stat<-numeric(newm)
    for (r in 1:newm){
      tot<-f[newR[r]+1]+3*log(300)+newneglog[r]
      stat[r]<-tot

    }
    f[t+1]<-min(stat)
    t1<-newR[which.min(stat)]
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

bestnormalize<-function(data,changepoints){
  datanew<-data
  changes<-length(changepoints)
  tau<-c(0,changepoints,length(data))
  for (c in 1:(changes+1)){
    datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]-mean(datanew[(tau[c]+1):tau[c+1]])
    subset<-datanew[(tau[c]+1):tau[c+1]]
    auto<-autocorrelations(subset)
    auto.iid<-whiteNoiseTest(auto,h0="iid",x=epsilon)
    pvalue<-auto.iid[["test"]][3]
    size<-length(subset)
    if (pvalue<0.05){
      dfw<-dfwhitenoise(subset)
      spline<-smooth.spline(subset,df=dfw)
      datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]-spline$y

    }
    datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]/sd(datanew[(tau[c]+1):tau[c+1]])

  }
  return(datanew)
}

###########
# Test for residualization and residualize
residualizebis<-function(met.data){
  resid<-FALSE
  plates <-
    unique(unlist(lapply(strsplit(
      rownames(met.data), "_"
    ), function(x) {
      return(x[2])
    })))
  if (length(plates)>=2){
    homo<-homogen.var(met.data)


    if (!homo){


      met.plates.sd <- c()
      print(plates)
      for (j in plates) {
        print(paste("plate:", j))
        grep.plates <-
          grep(paste("plate", "_", j, "_", sep = ""), row.names(met.data))
        sd.plate    <-  sd(met.data[grep.plates,], na.rm = T)

        print(paste("SD PLATE:", sd.plate))

        if(!is.na(sd.plate) & sd.plate == 0) {
          sd.plate <- NA
        }

        met.plates.sd <-
          c(met.plates.sd, rep(sd.plate, length(grep.plates)))
      }

      # Normalize
      met.data[["norm."]] <- met.data[["orig"]] / met.plates.sd


      met.data[["group"]] <-
        as.factor(paste("plate", "_", unlist(lapply(strsplit(rownames(met.data), "_"),
                                                    function(x) {return(x[2])})), sep = ""))

      # ANOVA test
      a <- aov(as.formula(paste("norm.", "~ group", sep = "")), data = met.data)
      a.pval <- summary(a)[[1]][, "Pr(>F)"][1]
      if (a.pval<0.05){
        resid<-TRUE
      }
    }else{
      met.data[["group"]] <-
        as.factor(paste("plate", "_", unlist(lapply(strsplit(rownames(met.data), "_"),
                                                    function(x) {return(x[2])})), sep = ""))
      met.data[["norm."]] <- met.data[["orig"]]
      # ANOVA test
      a <- aov(as.formula(paste("orig", "~ group", sep = "")), data = met.data)
      a.pval <- summary(a)[[1]][, "Pr(>F)"][1]
      if (a.pval<0.05){
        resid<-TRUE
      }
    }

    if (resid==TRUE){

      residualized<-residualize(met.data,"norm.","plate",plates)
      if (residualized$status == 0) {
        met.data[["resid"]]<-residualized$res
        return(met.data[["resid"]])

      }

    }else{
      return(met.data[["norm."]])
    }
  }
  return(met.data[["orig"]])
}
###########
# Normalization with finding change points and modified WINN

PELTwfknormalize<-function(data){
  corrected<-data
  for (k in 1:dim(data)[2]) {

  knots<-seq(1,length(data[[k]]),length.out=floor(length(data[[k]])/60)+2)
  knots<-knots[-c(1,length(knots))]
  knots<-knots+rep(0.5,length(knots))*(floor(knots)==knots)
  changepoints<-fkPELT(data[[k]],knots = knots)
  datanew<-changetometa(data[[k]],changepoints)
  changes<-length(changepoints)
  tau<-c(0,changepoints,length(data))
  datanew<-residualizebis(datanew)


  for (c in 1:(changes+1)){
    subset<-datanew[(tau[c]+1):tau[c+1]]
    auto<-autocorrelations(subset)
    auto.iid<-whiteNoiseTest(auto,h0="iid",x=epsilon,nlags=10,method="LjungBox")
    pvalue<-auto.iid[["test"]][3]
    size<-length(subset)
    if (!is.na(pvalue) & pvalue<0.01){
      if (size>=5){
        dfw<-dfwhitenoise(subset)
        spline<-smooth.spline(subset,df=dfw)
        datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]-spline$y
      }

    }


  }

  datanew<-changetometa(datanew,changepoints)

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

        p.wn <- Box.test(datanew[["orig"]], lag=20, type = "Ljung-Box")$p.value
        pass.wn<- p.wn>0.05
        corrected[,k]<-datanew[["orig"]]

      }

      p.wn <- Box.test(datanew[["norm."]], lag=20, type = "Ljung-Box")$p.value
      pass.wn<- p.wn>0.05

      corrected[,k]<-datanew[["norm."]]
    }
  }
  p.wn <- Box.test(datanew[["orig"]], lag=20, type = "Ljung-Box")$p.value
  pass.wn<- p.wn>0.05
  corrected[,k]<-datanew[["orig"]]
  }
  return(corrected)
}
########
# Compare method on real data
compare_real_data<-function(data){
  winn.pass.wn.2<-vector(length=dim(data)[2])
  winn.ch.pass.wn.2<-vector(length=dim(data)[2])
  f<-dim(data)[2]
  for (i in 1:f) {
    winn1<-winn(data[i],3,save.pdfs = FALSE)
    winn.ch<-changomics(data[[i]],1,save.pdfs = FALSE)
    winn.pass.wn.2[i]<-winn1$summary.transf["pass.wn.2"][[1]]
    winn.ch.pass.wn.2[i]<-winn.ch$summary.transf["pass.wn.2"][[1]]
  }

  percent.pass.ch<-mean(winn.ch.pass.wn.2)
  percent.pass.winn<-mean(winn.pass.wn.2)
  return(c(percent.pass.winn,percent.pass.ch))
}


########
# Test if find correctly changepoints
test<-function(nchange=6){
  x<-1:200
  y<-1:150
  #data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(146,mean=4,sd=1),rnorm(500,mean = sin(0.2*seq(1,125,by=0.25)),sd=0.5),rnorm(300,mean=seq(0.01,3,by=0.01)))
  #data<-c(rnorm(400,mean=0.01*x,sd=0.5),rnorm(300,mean =2-0.0001*y^2 ,sd=1),rnorm(200,mean = 0.0003*x^2,sd=1),rnorm(300,mean=0.5,sd=1.1),rnorm(300,mean = -0.5,sd=1.5),rnorm(100,mean = 0.5,sd=0.67))
  #data<-c(rnorm(450,mean = 0,sd=1),rnorm(150,mean=0.5,sd=1.1),rnorm(300,mean = 0,sd=1.5),rnorm(100,mean = 0.5,sd=1))
  #data<-rnorm(1000)
  tic()
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
  toc()
}

