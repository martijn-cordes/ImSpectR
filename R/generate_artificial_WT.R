#' @title Generate Artificial Wildtype
#'
#' @description Generate an Gaussian shaped Artificial Wildtype peak pattern
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples generate_artifical_WT(7)
#'
#' @export generate_artifical_WT
generate_artifical_WT <- function(no.peaks, plot=F, three.bp.positions=NULL, percentage_second_curve_height = NULL)
{
  #Generate peak pattern
  no.peaks <- no.peaks

  if (is.null(three.bp.positions)) {
    three.bp.positions <- 36
  }

  total_peak_width <- three.bp.positions*no.peaks

  x <- seq(0,total_peak_width, by = 1)
  mean <- mean(x)
  y <- dnorm(x, mean = mean, sd = total_peak_width/5)

  if (is.null(percentage_second_curve_height)) {
    percentage_second_curve_height <- 10
  }

  if(percentage_second_curve_height < 3) {
    percentage_second_curve_height <- 3
  }

  y_second_curve <- dnorm(x, mean = mean, sd = total_peak_width/5)
  y_second_curve <- y_second_curve * ((max(y) * percentage_second_curve_height / 100) / max(y_second_curve))

  position_peak <- c()
  for (i in 1:no.peaks-1) {
    position_peak <- c(position_peak, (round(three.bp.positions/2))+(three.bp.positions*i))
  }

  max_peak_no <- which(y[x %in% position_peak] == max(y[x %in% position_peak]))

  x_positions_single_peaks <- c()
  heights_single_peaks <- c()
  complete_peak_pattern <- c()
  slope <- c()
  sds <- c()

  for(p in 1:length(position_peak)) {
    mean_single_peak <- position_peak[p]
    x_single_peak <- seq(mean_single_peak-round(three.bp.positions/2),mean_single_peak+round(three.bp.positions/2), by = 1)

    current_slope <- (y[mean(which(round(x) == position_peak[p])[1])] - y_second_curve[mean(which(round(x) == x_single_peak[1])[1])])/(position_peak[p]- x_single_peak[1])
    slope <- c(slope,current_slope)


    #####
    #SD of heightest peak = width of three bp's / 5
    #SD of individual peaks is dependant on the slope of the heightest peak with above SD
    #Slope = y2-y1 / x2-x1
    #y2 = peak height
    #y1 = 0
    #x2 = position of peak
    #x1 = position of peak - SD
    #so x2-x1 = initial SD
    #new SD = height of peak to create / slope of heighest peak
    #new SD slope is right but SD is to small, add log10 of difference of height of peak to heightest peak to compensate
    #####


    sd <- round(three.bp.positions / (10-(log2(ceil(percentage_second_curve_height) +1)) ))
    slope_heightest_peak <- (y[mean(which(round(x) == position_peak[max_peak_no][1])[1])] / sd)
    sd <- y[mean(which(round(x) == position_peak[p])[1])] / slope_heightest_peak
    sd <- sd + (log2((max(y)/ y[mean(which(round(x) == position_peak[p])[1])]))*(log10(ceil(percentage_second_curve_height) +1)))
    sds <- c(sds, sd)

    height_single_peak <- y[mean(which(round(x) == position_peak[p])[1])]

    y_single_peak <- dnorm(x_single_peak, mean = mean_single_peak, sd = sd)

    if (p == 1 | p== 7) {
      y_single_peak  <- y_single_peak  * (height_single_peak / max(y_single_peak )/1.05)
    } else {
      y_single_peak  <- y_single_peak  * height_single_peak / max(y_single_peak )
    }

    if (percentage_second_curve_height > 1) {
      difference_in_lines <- y_second_curve[which(x %in% x_single_peak)] - y_single_peak

      overlap_1 <- order(difference_in_lines[1:which(diff(which(difference_in_lines > 0)) > 5)])[1]
      overlap_2 <- which(difference_in_lines > 0)[which(diff(which(difference_in_lines > 0)) > 5)+1]

      x_single_peak <- x_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]
      y_single_peak <- y_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]
    }


    complete_peak_pattern <- c(complete_peak_pattern, y_single_peak)
    x_positions_single_peaks <- c(x_positions_single_peaks, x_single_peak)
    height_single_peak <- max(y_single_peak)
    heights_single_peaks <- c(heights_single_peaks, height_single_peak)

  }

  cp <- c()
  cp[1:length(x)] <- NA
  cp[x_positions_single_peaks] <- complete_peak_pattern
  vec <- which(is.na(cp))
  groups <- split(vec, cumsum(seq_along(vec) %in% (which(diff(vec)>1)+1)))

  cp[c(1,length(cp))] <- 0
  cp <- approx(cp, xout=seq_along(cp))$y

  for (g in 1:length(groups)) {

    if(g == 1) {
      data <- cp[(unlist(groups[g])[1]) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct") # 10% smoothing span

      cp[(unlist(groups[g])[1]) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)] <- predict(loessMod10)
    }

    else if(g == length(groups)) {
      data <- cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))])]
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct") # 10% smoothing span
      cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))])] <- predict(loessMod10)
    }

    else {
      data <- cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct") # 10% smoothing span
      cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]  <- predict(loessMod10)

    }
  }


  complete_peak_pattern <- cp

  span <- percentage_second_curve_height / 100
  if(span > 0.06) {
    span <- 0.06
  }
  df <- data.frame(x=1:length(complete_peak_pattern), y=complete_peak_pattern)
  loessMod <- loess(y ~ x, data=df, span=span, surface="direct") # 6% smoothing span
  # get smoothed output
  artifical_WT <- predict(loessMod)

  if(plot == T) {
    plot(artifical_WT,type="l")
  }

  artifical_WT

}

