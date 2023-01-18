
#####
#' @title changomics.2
#' @description A function for metabolite correction. This function first estimate the change points
#'  (jumps) and then correct each different segment by detrending and normalizing when needed.
#' @param data input dataset of metabolites as data frame.
#' @param graph TRUE or FALSE. If TRUE will display the graph of uncorrected data with estimated change points and correction

#' @return the corrected dataframe of metabolite(s).
#' @import sarima
#' @import splines
#' @import dplyr
#' @importFrom stats na.omit
#' @export
#' @examples
#' met1<-rnorm(200)
#' changomics.2(as.data.frame(met1))

changomics.2<-function(data,graph=FALSE){
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
    #############
    #Detrending
    changes<-length(changepoints)
    datanew<-data[[k]]
    tau<-c(0,changepoints,length(data[[k]]))
    pred<-NULL
    for (c in 1:(changes+1)){
      subset<-datanew[(tau[c]+1):tau[c+1]]
      size<-length(subset)
      if (size>=5){
        pvalue<-Box.test(subset,lag=floor(length(subset)/2),type="Ljung-Box")$p.value
        if (!is.na(pvalue) & pvalue<0.05){
        dfw<-dfwhitenoise(subset)
        spline<-smooth.spline(subset,df=dfw)
        pred<-c(pred,spline$y)
      }else{
        pred<-c(pred,rep(mean(subset),size))
      }

      }else{
        pred<-c(pred,rep(mean(subset),size))
      }
    }
    if (k==1){
    if (graph==TRUE){
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
    datanew<-changetometa.2(as.data.frame(datanew),changepoints)
    colnames(datanew)<-"orig"
    if (changes>=1){

      homo2<-homogen.var.2(datanew)
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
          sd.plate    <-  sd(datanew[grep.plates,], na.rm = T)



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
  ret<-as.data.frame(data[[1]])
  colnames(ret)<-colnam
  row.names(ret)<-names
  return(ret)
}
# normalize<-function(data,changepoints){
#   changes<-length(changepoints)
#   datanew<-data
#   tau<-c(0,changepoints,length(data))
#   for (c in 1:(changes+1)){
#     subset<-datanew[(tau[c]+1):tau[c+1]]
#     size<-length(subset)
#
#     if (size>=5){
#       dfw<-dfwhitenoise(subset)
#       spline<-smooth.spline(subset,df=dfw)
#       datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]-spline$y
#     }else{
#       datanew[(tau[c]+1):tau[c+1]]<-datanew[(tau[c]+1):tau[c+1]]-mean(datanew[(tau[c]+1):tau[c+1]])
#     }
#
#
#
#   }
#   datanew<-changetometa.2(as.data.frame(datanew),changepoints)
#   colnames(datanew)<-"orig"
#   if (changes>=1){
#
#     homo2<-homogen.var.2(datanew)
#     if (!homo2){
#
#       plates <-
#         unique(unlist(lapply(strsplit(
#           rownames(datanew), "_"
#         ), function(x) {
#           return(x[2])
#         })))
#
#       met.plates.sd <- c()
#       print(plates)
#       for (j in plates) {
#         print(paste("plate:", j))
#         grep.plates <-
#           grep(paste("plate", "_", j, "_", sep = ""), row.names(datanew))
#         sd.plate    <-  sd(datanew[grep.plates,], na.rm = T)
#
#         print(paste("SD PLATE:", sd.plate))
#
#         if(!is.na(sd.plate) & sd.plate == 0) {
#           sd.plate <- NA
#         }
#
#         met.plates.sd <-
#           c(met.plates.sd, rep(sd.plate, length(grep.plates)))
#       }
#
#       # Normalize
#       datanew[["norm"]] <- datanew[["orig"]] / met.plates.sd
#       if (sum(is.na(datanew[["norm"]]))>0){
#
#         return(datanew[["orig"]])
#
#       }
#
#       return(datanew[["norm"]])
#     }
#   }
#
#   return(datanew[["orig"]])
#
# }

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
    value[i-1]<-Box.test(epsilon,lag=15,type="Ljung-Box")$p.value
  }
  df<-which.max(value)+1
  return(df)
}

homogen.var.2 <- function(met.dat)
{
  # Create a vector to store results of homogeneity variance test
  is.hom.var <- FALSE

    met.df <-  met.dat["orig"]

    if(sd(met.dat[,1], na.rm=T) == 0){
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
          fligner.test(as.formula(paste("orig", " ~ group", sep = "")), data = met.df)
        p.val <- flig.tst$p.value
      } else{
        # Change this once "car" library is available
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










