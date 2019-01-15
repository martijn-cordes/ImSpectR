#' @title Score in-silico pattern
#'
#' @description Score in-silico pattern, should get an score of 100. Used for generating score of all peak patterns
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples score_perfect_plot(three.bp.positions, no.peaks,percentage_second_curve_height)
#'
#' @export score_perfect_plot
score_perfect_plot <- function(three.bp.positions, no.peaks, percentage_second_curve_height, plot=F) {

  percentage = percentage_second_curve_height
  ########################################
  ########################################
  ####
  #perfect plot benchmark scoring

  sample_name <- "In Silico"
  #sample_index <- which(names(all_basepair_positions) == sample_name)
  total_peak_width <- three.bp.positions*no.peaks

  artifical_WT <- generate_artifical_WT(no.peaks, three.bp.positions = 36, percentage_second_curve_height=percentage, plot=F)
  artifical_WT <- artifical_WT*1000000

  pattern <-artifical_WT
  sample <- artifical_WT


  alignment <- dtw(pattern,sample, keep=TRUE,step=asymmetricP05, open.end=T,open.begin=T,window.type='sakoechiba',window.size=1)

  d <- alignment
  xts <- d$query
  idx <- 1:length(d$index1)

  if (plot == T) {
    plot(d$reference, type="l", xlim = c(d$index2[idx][1]-25 , d$index2[idx][length(idx)]+25), main="In silico control", xlab="", ylab="", yaxt="n", xaxt="n")
    #points(d$query, col="red", type="l")
    #points(d$index2[idx], xts[d$index1[idx]],type="l", col="darkgreen", lwd=2)
  }


  #matched <- xts[d$index1[idx]][-which(duplicated(d$index2[idx]))]
  matched <- xts[d$index1[idx]]


  #df <- data.frame(x=d$index2[idx][-which(duplicated(d$index2[idx]))], y=matched)
  df <- data.frame(x=d$index2[idx], y=matched)
  loessMod5 <- loess(y ~ x, data=df, span=0.05, surface="direct") # 5% smoothing span
  # get smoothed output
  smoothed5 <- predict(loessMod5)

  df_smooth <-df
  df_smooth$y <- smoothed5
  df_smooth$distance <- d$normalizedDistance

  alignment <- df_smooth

  peaks <- fpeaks(alignment$y, f=200, amp=c(200,200), plot=F)

  bp_positions_dtw_peaks <- c()
  x_positions_dtw_peaks <- c()
  peak_heights_dtw_peaks <- c()
  peak.margin <- 20

  for (p in 1:nrow(peaks)) {

    if (alignment$x[which(alignment$y %in% peaks[,2][p])] - peak.margin < 0) {
      peak.margin <- alignment$x[which(alignment$y %in% peaks[,2][p])]
    }

    if (alignment$x[which(alignment$y %in% peaks[,2][p])] + peak.margin > length(sample)) {
      peak.margin <- length(sample) - alignment$x[which(alignment$y %in% peaks[,2][p])]
    }

    if (alignment$x[which(alignment$y %in% peaks[,2][p])][1] %in% x_positions_dtw_peaks) {
      position <- alignment$x[which(alignment$y %in% peaks[,2][p])][2]
    } else {
      position <- alignment$x[which(alignment$y %in% peaks[,2][p])][1]
    }

    x_pos <- which(pattern %in% sample[which(sample == max(sample[(position-peak.margin) : (position+peak.margin)], na.rm=T))])
    x_positions_dtw_peaks <- c(x_positions_dtw_peaks,x_pos)

    peak_heights_dtw_peaks <- c(peak_heights_dtw_peaks,
                                sample[which(sample == max(sample[(position-peak.margin) : (position+peak.margin)], na.rm=T))]
    )

    bp_positions_dtw_peaks <- c(bp_positions_dtw_peaks, p*3)

  }

  if (plot == T) {
    points(x_positions_dtw_peaks, peak_heights_dtw_peaks, col="red")
  }

  score_curve_start <- 1
  score_curve_end <- length(sample)

  #plot(all_basepair_positions[[1]]$yy[pos_template], type="l")

  x <- seq(score_curve_start,score_curve_end, by = 1)
  mean <- mean(x)-1
  y <- dnorm(x, mean = mean, sd = total_peak_width/5.2)
  y <- y * (max(sample) / max(y))

  if (plot == T) {
    points(x,y, type="l")
  }

  ##SCORING
  #score <- as.numeric((trapz(sample) / trapz(artifical_WT))*100)

  position_peak <- x_positions_dtw_peaks
  max_peak_no <- which(y[x %in% position_peak] == max(y[x %in% position_peak]))


  x_pos_max_peak <- x_positions_dtw_peaks[which(peak_heights_dtw_peaks == max(peak_heights_dtw_peaks))]
  second_curve_height <- min(sample[(x_pos_max_peak-three.bp.positions) : (x_pos_max_peak+three.bp.positions)]) / max(peak_heights_dtw_peaks)

  y_second_curve <- dnorm(x, mean = mean, sd = total_peak_width/5)
  y_second_curve <- y_second_curve * (((max(y) * second_curve_height) / max(y_second_curve)) *1.1)

  if (plot == T) {
    points(x,y_second_curve, type="l", col="red")
  }

  heights_single_peaks <- c()
  heights_sample_peaks <- c()
  complete_peak_pattern <- c()
  slope <- c()
  sds <- c()
  peak_area_sample <- c()
  peak_area_wt <- c()
  peak_area_second_curve_single_peaks <- c()
  single_peak_penalties <- c()

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
    percentage_second_curve_height <- max(y_second_curve) / max(y) * 100

    sd <- round(three.bp.positions / (9-(log2(ceil(percentage_second_curve_height) +1)) ))
    slope_heightest_peak <- (y[mean(which(round(x) == position_peak[max_peak_no][1])[1])] / sd)
    sd <- y[mean(which(round(x) == position_peak[p])[1])] / slope_heightest_peak
    sd <- sd + (log2((max(y)/ y[mean(which(round(x) == position_peak[p])[1])]))*(log10(ceil(percentage_second_curve_height) +1)))
    sds <- c(sds, sd)

    height_single_peak <- y[mean(which(round(x) == position_peak[p])[1])]

    y_single_peak <- dnorm(x_single_peak, mean = mean_single_peak, sd = sd)
    y_single_peak  <- y_single_peak  * (height_single_peak -  (y_second_curve[mean(which(x %in% x_single_peak))])/2)     / max(y_single_peak )

    y_single_peak <- y_single_peak + (y_second_curve[mean(which(x %in% x_single_peak))]/2)

    if (percentage_second_curve_height > 3) {
      difference_in_lines <- y_second_curve[which(x %in% x_single_peak)] - y_single_peak

      overlap_1 <- order(difference_in_lines[1:which(diff(which(difference_in_lines > 0)) > 5)])[1]
      overlap_2 <- which(difference_in_lines > 0)[which(diff(which(difference_in_lines > 0)) > 5)+1]


      #overlaps <- order(abs(y_second_curve[which(x %in% x_single_peak)] - y_single_peak))
      #overlaps <- overlaps[c(1,(which(abs(diff(overlaps)) > 5) + 1)[1])]

      # overlaps <- overlaps[order(overlaps)]
      #overlap_1 <- overlaps[1]
      #overlap_2 <- overlaps[2]

      x_single_peak <- x_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]
      y_single_peak <- y_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]
    }


    complete_peak_pattern <- c(complete_peak_pattern, y_single_peak)
    #x_positions_single_peaks <- c(x_positions_single_peaks, x_single_peak)
    height_single_peak <- max(y_single_peak)
    heights_single_peaks <- c(heights_single_peaks, height_single_peak)
    height_sample_peak <-  max(sample[x_single_peak])
    heights_sample_peaks <- c(heights_sample_peaks, height_sample_peak)


    if(plot==T) {
      points(x_single_peak,y_single_peak, type="l", col="blue")
    }

    y_single_peak <- y_single_peak - y_second_curve[which(x %in% x_single_peak)]
    y_single_peak[which(y_single_peak < 0)] <- 0

    #points(x_single_peak, y_second_curve[which(x %in% x_single_peak)])

    corrected_sample <- sample[x_single_peak] - y_second_curve[which(x %in% x_single_peak)]
    corrected_sample[which(corrected_sample < 0)] <- 0


    #data outside of wt peak is penalized
    if (height_single_peak > height_sample_peak) {
      single_peak_penalty_area <- sum((corrected_sample - y_single_peak)[which((corrected_sample - y_single_peak) > 0 )])
    } else {
      single_peak_penalty_area <- 0
    }


    single_peak_penalties <- c(single_peak_penalties, single_peak_penalty_area)
    peak_area_corrected_sample <- as.numeric(trapz(corrected_sample))
    #peak_area_corrected_sample <- as.numeric(trapz(corrected_sample)) - single_peak_penalty_area

    peak_area_wt <- c(peak_area_wt, trapz(y_single_peak))
    peak_area_sample <- c(peak_area_sample, peak_area_corrected_sample)
    peak_area_second_curve_single_peaks <- c(peak_area_second_curve_single_peaks, trapz(y_second_curve[which(x %in% x_single_peak)]))



  }
  #######
  #calculate percentage difference between the peaks
  fraction <- peak_area_wt / (peak_area_sample)
  if(is.infinite(max(fraction))) {
    fraction[which(is.infinite(fraction))] <- sds[is.infinite(fraction)]
  }

  #wt peak area smaller then measured data area
  fraction[which(fraction < 1)] <- peak_area_wt[which(fraction < 1)] / peak_area_sample[which(fraction < 1)]
  #less penalty for peaks bigger then curve because this adds complexity
  fraction[which(fraction < 1)] <- (peak_area_sample[which(fraction < 1)] / peak_area_wt[which(fraction < 1)])

  plot_area <-trapz(sample)
  bell_curve_area <- trapz(y)
  second_curve_area <- trapz(y_second_curve)


  #0.9360127
  within_correction <- (plot_area - second_curve_area) / sum(peak_area_sample)
  #within_correction <- 0.9360127

  area_in_between_peaks <- (plot_area - second_curve_area) - (sum(peak_area_sample) * within_correction)
  #area_in_between_peaks <- area_in_between_peaks + sum(single_peak_penalties)
  max_area_in_between_wt_peaks <- (bell_curve_area - second_curve_area) - sum(peak_area_wt)

  percentage_in_between <- area_in_between_peaks / (max_area_in_between_wt_peaks / 100)

  plot_score <- (length(fraction) / sum(fraction)) * 100
  plot_score <- plot_score - (plot_score*((percentage_in_between*2)/100))
  plot_score_correction <- 100/plot_score
  #plot_score_correction <- 1.037922
  plot_score <- plot_score * plot_score_correction

  if (plot == T) {
    legend("topright",
           paste("Peak score :", round(plot_score,2))
    )
  }


  score_adjustements <- c(within_correction, plot_score_correction)

}
