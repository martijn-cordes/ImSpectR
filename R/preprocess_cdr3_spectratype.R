#' @title Preprocess spectratype data
#'
#' @description Preprocess spectratype data into an score-ready object
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples preprocess_cdr3_spectratype(folder)
#'
#' @export get_basepair_positions
preprocess_cdr3_spectratype <- function (folder)
{
  
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  
  if (substr(folder, nchar(folder),nchar(folder)) != "/") {
    folder <- paste(folder,"/",sep="")
  }
  
  my.ladder <- c(35, 50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400, 450, 490, 500)
  fsa.loaded <- storing.inds(folder)
  spectratype_dataset <- get_basepair_positions(fsa.loaded, cols = 1, my.ladder, channel.ladder=NULL,  init.thresh=1750, ladd.init.thresh=1000)
  
  spectratype_dataset
  
}