#' @title Align peak pattern
#'
#' @description Perform Dynamic Time warping to match an in-silico peak pattern to user data
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples align_peak_pattern(in_silico_peak_pattern, sample_data)
#'
#' @export align_peak_pattern
align_peak_pattern <- function(pattern, sample, plot.pattern.matching=F, start=NULL, pattern.width=NULL, window.size=1, sliding.window=T)
{
  #pattern <- artifical_WT
  #sample <- generate_monoclonal_artifical_WT(10, fsa.names[1])

  if(is.null(sliding.window)) {
    sliding.window <- T
  }

  all_basepair_positions <- sample
  #i <- 1

  if (sliding.window==T) {
    window_start <- which(all_basepair_positions[[1]]$xx > start)[1]
    window_end <- which(all_basepair_positions[[1]]$xx > start+pattern.width)[1]
    sample <- all_basepair_positions[[1]]$yy[window_start:window_end]
  }


  if (sliding.window==T) {
    alignment <- dtw(pattern,sample, keep=T,step=asymmetricP05, open.end=T,open.begin=F,window.type='sakoechiba',window.size=window.size)
  } else {
    alignment <- dtw(pattern,sample, keep=TRUE,step=asymmetricP05, open.end=T,open.begin=T)
  }

  #alignment <- dtw(pattern,sample, keep=TRUE,step=asymmetricP2, open.end=T,open.begin=T,window.type='sakoechiba',window.size=1)
  #alignment <- dtw(pattern,sample, keep=TRUE,step=symmetricP05, window.size=1)

  d <- alignment
  xts <- d$query
  idx <- 1:length(d$index1)
  #plot(d$reference, type="l", xlim = c(d$index2[idx][1]-100 , d$index2[idx][length(idx)]+100), main=fsa.names[i])
  #points(d$query, col="red", type="l")
  #points(d$index2[idx], xts[d$index1[idx]],type="l", col="darkgreen", lwd=2)

  matched <- xts[d$index1[idx]][-which(duplicated(d$index2[idx]))]


  df <- data.frame(x=d$index2[idx][-which(duplicated(d$index2[idx]))], y=matched)
  loessMod5 <- loess(y ~ x, data=df, span=0.05, surface="direct") # 5% smoothing span
  # get smoothed output
  smoothed5 <- predict(loessMod5)
  pattern_range <- df$x + window_start

  df_smooth <-df
  df_smooth$y <- smoothed5
  df_smooth$distance <- d$normalizedDistance
  df_smooth$bp <- all_basepair_positions[[1]]$xx[pattern_range]
  df_smooth$x_pos <- pattern_range



  if(plot.pattern.matching == T) {
    if(sliding.window == T) {
      plot(all_basepair_positions[[1]]$yy , type="l", main=names(all_basepair_positions)[1], xlim=c(window_start-100, window_end+100),ylim=c(0,max(sample)), xaxt="n")
    } else {
      plot(all_basepair_positions[[1]]$yy , type="l", main=names(all_basepair_positions)[1], xlim=c(df_smooth$x[1], df_smooth$x[length(df_smooth$x)]),ylim=c(0,max(df_smooth$y)), xaxt="n")
    }

    ticks <- seq(pattern_range[1]-100, pattern_range[length(pattern_range)]+100, by = 20)
    labels <- round(all_basepair_positions[[1]]$xx[ticks])

    axis(side=1, at=ticks, labels=labels)

    if(sliding.window==T) {
      points(df_smooth$x + window_start,df_smooth$y,type="l", col="darkgreen", lwd=2)
    } else {
      points(df_smooth$x,df_smooth$y,type="l", col="darkgreen", lwd=2)
    }

    legend("topright", paste("Distance score: ", round(df_smooth$distance[1],2)))

  }

  df_smooth
}


