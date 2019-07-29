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
  
  sample_name <- "In Silico"
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
  }
  
  matched <- xts[d$index1[idx]]
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
  
  x <- seq(score_curve_start,score_curve_end, by = 1)
  mean <- mean(x)-1
  y <- dnorm(x, mean = mean, sd = total_peak_width/5.05)
  y <- y * (max(sample) / max(y))
  
  if (plot == T) {
    points(x,y, type="l")
  }
  
  position_peak <- x_positions_dtw_peaks
  max_peak_no <- which(y[x %in% position_peak] == max(y[x %in% position_peak]))
  
  
  #Fit  all normal curves
  min_heights <- c()
  min_positions <- c()
  for(pp in 1:length(position_peak)) {
    if(pp == 1) {
      min_position <- which(x %in% (position_peak[pp]-10)) : which(x %in% (position_peak[pp]+(three.bp.positions+10)))
    }
    else if(pp == length(position_peak)) {
      min_position <- which(x %in% (position_peak[pp]-10)) : which(x %in% (position_peak[pp]+10)) 
    } else {
      min_position <- which(x %in% (position_peak[pp]-(three.bp.positions-10))) : which(x %in% (position_peak[pp]+(three.bp.positions+10)))
    }
    
    min_positions <- c(min_positions, list(min_position))
    
    minimum <- min(sample[min_position])
    min_heights <- c(min_heights, minimum)
  }
  
  second_curve_height <- min_heights[round(length(min_heights)/2)] / max(peak_heights_dtw_peaks)
  
  
  x_pos_max_peak <- x_positions_dtw_peaks[which(peak_heights_dtw_peaks == max(peak_heights_dtw_peaks))]
  second_curve_height <- min(sample[(x_pos_max_peak-three.bp.positions) : (x_pos_max_peak+three.bp.positions)]) / max(peak_heights_dtw_peaks)
  
  y_second_curve <- dnorm(x, mean = mean, sd = total_peak_width/6)
  y_second_curve <- y_second_curve * (((max(y) * second_curve_height) / max(y_second_curve)) *1)
  
  if((min_heights[2] < min(y_second_curve[min_positions[[2]] ]))) {
    y_second_curve <- y_second_curve / (min(y_second_curve[min_positions[[2]]]) / min_heights[2])
  }
  
  if ((max(y_second_curve) / max(y)) > 0.3 | (min_heights[2] < min(y_second_curve[min_positions[[2]] ]))) {
    y_second_curve <- y_second_curve * (((max(y) * 0.05) / max(y_second_curve)) *1.5)
  }
  
  
  y_second_curve <- y_second_curve-min(y_second_curve)
  
  y_second_curve <- dnorm(x, mean = mean, sd = total_peak_width/6)
  y_second_curve <- y_second_curve * (((max(y) * second_curve_height) / max(y_second_curve)) *1)
  
  
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
  x_positions_single_peaks <-c()
  all_overlaps <- c()
  
  for(p in 1:length(position_peak)) {
    mean_single_peak <- position_peak[p]
    x_single_peak <- seq(mean_single_peak-round(three.bp.positions/2),mean_single_peak+round(three.bp.positions/2), by = 1)
    
    current_slope <- (y[mean(which(round(x) == position_peak[p])[1])] - y_second_curve[mean(which(round(x) == x_single_peak[1])[1])])/(position_peak[p]- x_single_peak[1])
    slope <- c(slope,current_slope)
    
    percentage_second_curve_height <- max(y_second_curve) / max(y) * 100
    
    sd <- round(three.bp.positions / (9-(log2(ceil(percentage_second_curve_height) +1)) ))
    slope_heightest_peak <- (y[mean(which(round(x) == position_peak[max_peak_no][1])[1])] / sd)
    sd <- y[mean(which(round(x) == position_peak[p])[1])] / slope_heightest_peak
    sd <- sd + (log2((max(y)/ y[mean(which(round(x) == position_peak[p])[1])]))*(log10(ceil(percentage_second_curve_height) +1)))
    sds <- c(sds, sd)
    
    height_single_peak <- y[mean(which(round(x) == position_peak[p])[1])]
    
    y_single_peak <- dnorm(x_single_peak, mean = mean_single_peak, sd = sd)
    y_single_peak  <- y_single_peak  * height_single_peak / max(y_single_peak )
    
    overlaps <- order(abs(y_second_curve[which(x %in% x_single_peak)] - y_single_peak))
    overlaps <- overlaps[c(1,(which(abs(diff(overlaps)) > 7) + 1)[1])]
    
    overlaps <- overlaps[order(overlaps)]
    overlap_1 <- overlaps[1]
    overlap_2 <- overlaps[2]
    
    
    all_overlaps <- c(all_overlaps,list(x_single_peak[overlaps]))

    #if ((max(y_second_curve) / max(y) * 100) > 10 ) {
    
    x_single_peak <- x_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]
    y_single_peak <- y_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]
    
    #}
    
    complete_peak_pattern <- c(complete_peak_pattern, y_single_peak)
    
    
    x_positions_single_peaks <- c(x_positions_single_peaks, x_single_peak)
    
    height_single_peak <- max(y_single_peak)
    heights_single_peaks <- c(heights_single_peaks, height_single_peak)
    height_sample_peak <-  max(sample[x_single_peak])
    heights_sample_peaks <- c(heights_sample_peaks, height_sample_peak)
    
    
    if(plot==T) {
      points(x_single_peak,y_single_peak, type="l", col="blue")
    }
    
    
    y_single_peak <- y_single_peak - y_second_curve[which(x %in% x_single_peak)]
    y_single_peak[which(y_single_peak < 0)] <- 0
    
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
    peak_area_wt <- c(peak_area_wt, trapz(y_single_peak))
    peak_area_sample <- c(peak_area_sample, peak_area_corrected_sample)
    peak_area_second_curve_single_peaks <- c(peak_area_second_curve_single_peaks, trapz(y_second_curve[which(x %in% x_single_peak)]))
    
  }
  
  #######
  #calculate percentage difference between the peaks
  peak_area_sample[which(peak_area_sample < (max(peak_area_sample) * 0.01))] <- max(peak_area_sample) * 0.01
  
  
  fraction <- peak_area_wt / (peak_area_sample)
  if(is.infinite(max(fraction))) {
    fraction[which(is.infinite(fraction))] <- sds[is.infinite(fraction)]
  }
  
  #wt peak area smaller then measured data area
  fraction[which(fraction < 1)] <- peak_area_wt[which(fraction < 1)] / peak_area_sample[which(fraction < 1)]
  #less penalty for peaks bigger then curve because this adds complexity
  fraction[which(fraction < 1)] <- (peak_area_sample[which(fraction < 1)] / peak_area_wt[which(fraction < 1)])
  
  
  #Fraction of difference of peak heights vs bell curve
  peak_heights_dtw_peaks[which(peak_heights_dtw_peaks < (max(peak_heights_dtw_peaks) * 0.01))] <- max(peak_heights_dtw_peaks) * 0.01
  
  fraction_heights <- peak_heights_dtw_peaks / heights_single_peaks
  #wt peak area smaller then measured data area
  fraction_heights[which(fraction_heights < 1)] <- peak_heights_dtw_peaks[which(fraction_heights < 1)] / heights_single_peaks[which(fraction_heights < 1)]
  #less penalty for peaks bigger then curve because this adds complexity
  fraction_heights[which(fraction_heights < 1)] <- (heights_single_peaks[which(fraction_heights < 1)] / peak_heights_dtw_peaks[which(fraction_heights < 1)])
  
  if(is.infinite(max(fraction_heights))) {
    fraction_heights[which(is.infinite(fraction_heights))] <- 100
  }
  
  
  
  
  ###Baseline correction 
  sample[x] <- sample[x] - y_second_curve
  sample[x][which(sample[x] < 0)] <- 0
  y <- y - y_second_curve

  in_between_seperate <- c()
  in_between_max <- c()
  
  in_between_seperate <- c(in_between_seperate, sum(sample[(all_overlaps[[1]][1]-5):all_overlaps[[1]][1]]))
  in_between_max <- c(in_between_max, sum(y[(which(x == all_overlaps[[1]][1])-5) : which(x == all_overlaps[[1]][1])]))
    

  for(o in 1:(length(all_overlaps)-1)) {
    in_between_seperate <- c(in_between_seperate, sum(sample[all_overlaps[[o]][2]:all_overlaps[[o+1]][1]]))
    in_between_max <- c(in_between_max, sum(y[which(x == all_overlaps[[o]][2]) : which(x == all_overlaps[[o+1]][1])]))
  }
  
  in_between_seperate <- c(in_between_seperate, sum(sample[all_overlaps[[length(all_overlaps)]][2]:(all_overlaps[[length(all_overlaps)]][2]+5)]))
  in_between_max <- c(in_between_max, sum(y[which(x == all_overlaps[[length(all_overlaps)]][2]) : (which(x == all_overlaps[[length(all_overlaps)]][2]+5))])) 
  
  
  #Score calculation
  plot_score <- (length(fraction) / sum(fraction)) * 100
  peak_score_correction <- 100/plot_score

  score_1 <- plot_score * peak_score_correction
  
  plot_score_heights <- (length(fraction_heights) / sum(fraction_heights)) * 100
  gaussian_correction <- 100 / plot_score_heights
  plot_score_heights <- plot_score_heights * gaussian_correction
  
  residual_penalty <- (trapz(in_between_seperate) / trapz(in_between_max)*100)
  residual_correction <- residual_penalty
  residual_penalty <- residual_penalty - residual_correction
  
  combined_score <- ((score_1+plot_score_heights)/2) -  residual_penalty

  
  if (plot == T) {
    legend("topright",
           paste("Peak score :", round(combined_score,2))
    )
  }
  
  score_adjustements <- c(peak_score_correction, gaussian_correction, residual_correction)
  
}
