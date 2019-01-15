#' @title Get basepair positions
#'
#' @description Calculate every basepair position for every datapoint in the loaded fsa files
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples get_basepair_positions(fsa.loaded, cols=1, ladder=ladder, init.thresh=200)
#'
#' @export get_basepair_positions
get_basepair_positions <- function (my.inds, cols = 1, ladder, channel.ladder = NULL, init.thresh = 200, ladd.init.thresh = 200)
{

  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)

  ci.upp = 1.96
  ci.low = 1.96
  dev = 50
  warn = FALSE
  window = 0.5
  method = "iter2"
  env = parent.frame()
  pref = 3


  if (is.null(channel.ladder)) {
    channel.ladder <- dim(my.inds[[1]])[2]
  }
  else {
    channel.ladder <- channel.ladder
  }
  if (dim(my.inds[[1]])[2] < channel.ladder) {
    print(paste("ERROR MY FRIEND!! you have indicated an argument channel.ladder=5, but your data contains less channels/colors"))
    stop
  }

  n.inds <- c(1:length(my.inds))

  count <- 0
  tot <- length(n.inds)
  pb <- txtProgressBar(style = 3)
  setTxtProgressBar(pb, 0)
  my.inds2 <- list(NA)
  thresh2 <- list(NA)
  for (i in 1:length(n.inds)) {
    count <- count + 1
    v1 <- n.inds[i]
    my.inds2[[i]] <- my.inds[[v1]]
    names(my.inds2)[i] <- names(my.inds)[i]
  }
  list.data <- list(NA)
  if (exists("list.data.covarrubias")) {
    list.data <- env$list.data.covarrubias
  }
  else {
    list.ladders <- lapply(my.inds2, function(x) {
      y <- x[, channel.ladder]
      return(y)
    })
    list.data <- lapply(list.ladders, find.ladder, ladder = ladder,
                        draw = F, dev = dev,
                        warn = warn, method = method, init.thresh = ladd.init.thresh)
  }
  list.models <- lapply(list.data, function(da) {
    y <- da[[3]]
    x <- da[[1]]
    mod <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5),
              data = da)
    return(mod)
  })
  list.models.inv <- lapply(list.data, function(da) {
    x <- da[[3]]
    y <- da[[1]]
    mod <- lm(y ~ x, data = da)
    return(mod)
  })
  xx <- lapply(my.inds2, function(x, cols) {
    1:length(x[, cols])
  }, cols = cols)
  newxx <- numeric()
  newyy <- numeric()
  new.whole.data <- list(NA)
  for (h in 1:length(xx)) {
    h1 <- n.inds[h]
    count <- count + 1
    newxx <- as.vector(predict(list.models[[h1]], newdata = data.frame(x = xx[[h]])))
    newyy <- my.inds2[[h]][, cols]
    new.whole.data[[h]] <- list(xx = newxx, yy = newyy)
    setTxtProgressBar(pb, (count/tot) * 0.25)
  }

  all_basepair_positions <- new.whole.data

  ##Only select linear part of bp conversion curve, cutoff data after shift of slope by calculating derivative of curve
  for(i in 1:length(all_basepair_positions)) {
    bp_data <- data.frame(x=1:length(all_basepair_positions[[i]]$xx), y=all_basepair_positions[[i]]$xx)
    #plot(y ~ x, data = as.data.frame(bp_data), main=names(all_samples)[i])

    DeltaY <- diff(bp_data$y)
    Turns <- which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) + 1
    #print(Turns)

    if(length(Turns)>1) {
      if (diff(Turns)>5000) {
        all_basepair_positions[[i]]$xx <- all_basepair_positions[[i]]$xx[Turns[1]:Turns[2]]
        all_basepair_positions[[i]]$yy <- all_basepair_positions[[i]]$yy[Turns[1]:Turns[2]]
      }
      #points(bp_data$x[Turns], bp_data$y[Turns], pch=16, col="red")
    }
    if(length(Turns)==1) {
      if ((Turns)>6000) {
        all_basepair_positions[[i]]$xx <- all_basepair_positions[[i]]$xx[0:Turns[1]]
        all_basepair_positions[[i]]$yy <- all_basepair_positions[[i]]$yy[0:Turns[1]]
      }
      #points(bp_data$x[Turns], bp_data$y[Turns], pch=16, col="red")
    }


    #plot(all_basepair_positions[[i]]$xx, main=names(all_samples)[i])
  }

  names(all_basepair_positions) <- gsub(".fsa", "",names(my.inds))

  all_basepair_positions

}

