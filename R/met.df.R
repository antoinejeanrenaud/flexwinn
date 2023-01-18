#' met.df
#'
#' @docType data
#'
#' @usage data(met.df)
#'
#' @format an object of class data.frame with 2 simluated metabolites data
#'
#'
#'\describe{ In this dataset we are giving two metabolites data that we
#'generated ourselves. "met1" is just some white noise from a normal
#'distribution and "met2" is a metaboite with jumps and trends to correct}
#'
#' @keywords dataset
#'
#' @examples
#' \dontrun{
#' data(met.df)
#' met.df<-as.data.frame(met.df)
#' corrected<-changomics(met.df)
#'}
#'
"met.df"
