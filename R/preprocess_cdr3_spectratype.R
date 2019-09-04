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

  #Preprocess
  fsa.loaded <- storing.inds(folder)
  my.ladder <- c(35, 50, 75, 100, 139, 150, 160, 200, 250,300, 340, 350, 400, 450, 490, 500)

  for (i in 1:length(fsa.loaded)) {

    fsa.loaded[[i]][,1][which(fsa.loaded[[i]][,1] < 0)] <- 0
    fsa.loaded[[i]][,2][which(fsa.loaded[[i]][,2] < 0)] <- 0
    fsa.loaded[[i]][,4][which(fsa.loaded[[i]][,4] < 0)] <- 0


    left <- which(fsa.loaded[[i]][,4] > (max(fsa.loaded[[i]][,4][3000:length(fsa.loaded[[i]][,4])])+100)) - 10
    right <- which(fsa.loaded[[i]][,4] > (max(fsa.loaded[[i]][,4][3000:length(fsa.loaded[[i]][,4])])+100)) + 10
    for(d in 1:length(left)) {

      fsa.loaded[[i]][,4][left[d]:right[d]] <- 0

    }

  }


  ladder.info.attach(stored=fsa.loaded, ladder=my.ladder, ladd.init.thresh=300, draw=F,method="iter2")
  spectratype_dataset <- get_basepair_positions(fsa.loaded, cols = 1, my.ladder, channel.ladder=NULL,  init.thresh=1750, ladd.init.thresh=300)

  spectratype_dataset

}
