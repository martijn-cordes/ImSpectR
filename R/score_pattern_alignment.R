#' @title Score peak pattern
#'
#' @description Score peak pattern
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples score_pattern_alignment(alignment, sample)
#'
#' @export score_pattern_alignment
score_pattern_alignment <- function(bp_positions_dtw_peaks, no.peaks, three.bp.positions, total_peak_width, plot.curve.fitting=T, plot=T, alt.scores=F)
{
  sample <- all_basepair_positions[[1]]$yy[pos_template]

  score_curve_start <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)[1]
  score_curve_start <- score_curve_start - round(three.bp.positions / 1)
  if(score_curve_start < 0) {
    score_curve_start <- 1
  }
  score_curve_end <- (score_curve_start+total_peak_width) + round(three.bp.positions / 1)

  if (all_basepair_positions[[1]]$xx[pos_template][score_curve_end] < bp_positions_dtw_peaks[length(bp_positions_dtw_peaks)]) {
    score_curve_end <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks[length(bp_positions_dtw_peaks)])+ round(three.bp.positions / 2)
  }

  x <- seq(score_curve_start-100,score_curve_end, by = 1)
  #assign("x", x, envir = .GlobalEnv)

  mean <- mean(x)
  mean <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)[4]

  sd <- total_peak_width/5.3
  #sd <- (which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)[4] -
  #  which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)[3]) * 2

  y <- dnorm(x, mean = mean, sd = sd)
  y <- y * (max(sample) / max(y))

  if(plot.curve.fitting==T) {  points(x,y, type="l") }

  position_peak <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)
  max_peak_no <- which(y[x %in% position_peak] == max(y[x %in% position_peak]))

  peak_heights_dtw_peaks <- c()
  for (p in 1:length(bp_positions_dtw_peaks)) {
    peak_heights_dtw_peaks <- c(peak_heights_dtw_peaks,
                                all_basepair_positions[[1]]$yy[pos_template][which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks[p] )]
    )
  }

  #Fit best normal curve
  #tab <- data.frame(position_peak=position_peak, height_peak=as.numeric(peak_heights_dtw_peaks))
  #mu <- tab$position_peak[which(tab$height_peak == max(tab$height_peak))]
  #k <-  tab$height_peak[which(tab$height_peak == max(tab$height_peak))]
  #Apply function nls
  #res <- nlsLM( height_peak ~ k*exp(-1/2*(position_peak-mu)^2/sigma^2), start=c(mu=mu,sigma=mu/10,k=k) , data = tab)
  #v <- summary(res)$parameters[,"Estimate"]
  #plot(function(position_peak) v[3]*exp(-1/2*(position_peak-v[1])^2/v[2]^2),col="red",add=T,xlim=c(position_peak[1]-200,position_peak[length(position_peak)]+200) )
  #bell_curve_points <- v[3]*exp(-1/2*(position_peak-v[1])^2/v[2]^2)
  #points(position_peak,bell_curve_points,col="red",add=T,xlim=c(tab$position_peak[1]-200,tab$position_peak[length(position_peak)]+200) )
  #x <- seq(score_curve_start,score_curve_end, by = 1)
  #mean <- v[1]
  #y_cor <- dnorm(x, mean = mean, sd = v[2])
  #y_cor <- y_cor * (v[3] / max(y_cor))
  #points(x,y_cor, type="l", col="blue")
  #y <- y_cor

  min_heights <- c()
  min_positions <- c()
  for(pp in 1:length(position_peak)) {
    min_position <- which(x %in% (position_peak[pp]-(three.bp.positions-10))) : which(x %in% (position_peak[pp]+(three.bp.positions+10)))
    min_positions <- c(min_positions, list(min_position))

    minimum <-min(all_basepair_positions[[1]]$yy[pos_template][(position_peak[pp]-(three.bp.positions-10)) : (position_peak[pp]+(three.bp.positions+10))])
    min_heights <- c(min_heights, minimum)
  }

  assign("min_heights", min_heights, envir = .GlobalEnv)

  second_curve_height <- min_heights[round(length(min_heights)/2)] / max(peak_heights_dtw_peaks)

  y_second_curve <- dnorm(x, mean = mean, sd = total_peak_width/4)
  y_second_curve <- y_second_curve * (((max(y) * second_curve_height) / max(y_second_curve)) *1.2)


  if((min_heights[2] < min(y_second_curve[min_positions[[2]] ]))) {
    y_second_curve <- y_second_curve / (min(y_second_curve[min_positions[[2]]]) / min_heights[2])

    #y_second_curve <- ((y_second_curve / (min(y_second_curve[min_positions[[2]]]) / min_heights[2])) +
    # (y_second_curve / (min(y_second_curve[min_positions[[6]]]) / min_heights[6]))  /2)
  }

  if ((max(y_second_curve) / max(y)) > 0.3 | (min_heights[2] < min(y_second_curve[min_positions[[2]] ]))) {
    y_second_curve <- y_second_curve * (((max(y) * 0.05) / max(y_second_curve)) *1.5)
  }


  y_second_curve <- y_second_curve-min(y_second_curve)
  assign("y_second_curve", y_second_curve, envir = .GlobalEnv)


  if(plot.curve.fitting == T) {
    points(x,y_second_curve, type="l", col="red")
  }


  #########################NEW SCORING


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

  for(p in 1:length(position_peak)) {
    mean_single_peak <- position_peak[p]

    x_single_peak <- seq(mean_single_peak-20,mean_single_peak+20, by = 1)
    #x_single_peak <- seq(mean_single_peak-round(three.bp.positions/2),mean_single_peak+round(three.bp.positions/2), by = 1)

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
    #assign("percentage_second_curve_height", percentage_second_curve_height, envir = .GlobalEnv)

    if(percentage_second_curve_height == 0) {
      percentage_second_curve_height <- 5
    }

    sd <- round(three.bp.positions / (10-(log2(ceil(percentage_second_curve_height) +1)) ))
    slope_heightest_peak <- (y[mean(which(round(x) == position_peak[max_peak_no][1])[1])] / sd)
    sd <- y[mean(which(round(x) == position_peak[p])[1])] / slope_heightest_peak
    sd <- sd + (log2((max(y)/ y[mean(which(round(x) == position_peak[p])[1])]))*(log10(ceil(percentage_second_curve_height) +1)))
    sds <- c(sds, sd)

    height_single_peak <- y[mean(which(round(x) == position_peak[p])[1])]

    y_single_peak <- dnorm(x_single_peak, mean = mean_single_peak, sd = sd)
    y_single_peak  <- y_single_peak  * height_single_peak / max(y_single_peak )

    if(plot==T) {
      #points(x_single_peak,y_single_peak, type="l", col="blue")
    }

    #if (percentage_second_curve_height > 1.5) {

    overlaps <- order(abs(y_second_curve[which(x %in% x_single_peak)] - y_single_peak))
    overlaps <- overlaps[c(1,(which(abs(diff(overlaps)) > 5) + 1)[1])]

    overlaps <- overlaps[order(overlaps)]
    overlap_1 <- overlaps[1]
    overlap_2 <- overlaps[2]

    x_single_peak <- x_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]

    #abline(v=x_single_peak[1])
    #abline(v=x_single_peak[length(x_single_peak)])

    y_single_peak <- y_single_peak[-c(c((1:overlap_1-1)),c((overlap_2+1):length(y_single_peak)))]

    if(length(duplicated(x_single_peak)) >0){
      # y_single_peak <- y_single_peak[-duplicated(x_single_peak)]
      #x_single_peak <- x_single_peak[-duplicated(x_single_peak)]

    }

    # }


    complete_peak_pattern <- c(complete_peak_pattern, y_single_peak)

    x_positions_single_peaks <- c(x_positions_single_peaks, x_single_peak)
    assign("x_positions_single_peaks", x_positions_single_peaks, envir = .GlobalEnv)

    height_single_peak <- max(y_single_peak)
    heights_single_peaks <- c(heights_single_peaks, height_single_peak)
    height_sample_peak <-  max(all_basepair_positions[[1]]$yy[pos_template][x_single_peak])
    heights_sample_peaks <- c(heights_sample_peaks, height_sample_peak)


    if(plot.curve.fitting==T) {
      points(x_single_peak,y_single_peak, type="l", col="blue")
    }

    #y_single_peak <- y_single_peak
    y_single_peak <- y_single_peak - y_second_curve[which(x %in% x_single_peak)]
    y_single_peak[which(y_single_peak < 0)] <- 0

    #points(x_single_peak, y_second_curve[which(x %in% x_single_peak)])

    #corrected_sample <- all_basepair_positions[[1]]$yy[pos_template][x_single_peak]
    corrected_sample <- all_basepair_positions[[1]]$yy[pos_template][x_single_peak] - y_second_curve[which(x %in% x_single_peak)]
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
    #peak_area_second_curve_single_peaks <- peak_area_second_curve_single_peaks*0


  }


  #####
  #####
  ####Plot overlay
  complete_peak_pattern <- c(rep(0,10),complete_peak_pattern, rep(0,10))

  x_positions_single_peaks <-
    c((x_positions_single_peaks[1]-10):(x_positions_single_peaks[1]),
      x_positions_single_peaks,
      (x_positions_single_peaks[length(x_positions_single_peaks)]):
        (x_positions_single_peaks[length(x_positions_single_peaks)]+10)
    )

  x_overlay<- 0:(x_positions_single_peaks[length(x_positions_single_peaks)] - x_positions_single_peaks[1])
  x_overlay_positions<- seq(x_positions_single_peaks[1],x_positions_single_peaks[length(x_positions_single_peaks)])

  x_positions_single_peaks_real <- x_positions_single_peaks - x_positions_single_peaks[1]

  cp <- c()
  cp[1:length(x_overlay)] <- NA
  cp[x_positions_single_peaks_real] <- complete_peak_pattern
  vec <- which(is.na(cp))
  groups <- split(vec, cumsum(seq_along(vec) %in% (which(diff(vec)>1)+1)))


  cp[c(1,length(cp))] <- 0
  cp <- approx(cp, xout=seq_along(cp))$y

  #assign("cp", cp, envir = .GlobalEnv)

  for (g in 1:length(groups)) {

    if(g == 1) {
      data <- cp[(unlist(groups[g])[1]) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]
      #plot(data)
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct") # 100% smoothing span
      #points(predict(loessMod10), type="l")
      cp[(unlist(groups[g])[1]) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)] <- predict(loessMod10)
    }

    else if(g == length(groups)) {
      data <- cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))])]

      #plot(data)
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct") # 100% smoothing span
      #points(predict(loessMod10), type="l")
      cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))])] <- predict(loessMod10)
    }

    else {
      data <- cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]

      #plot(data)
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct") # 100% smoothing span
      #points(predict(loessMod10), type="l")
      cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]  <- predict(loessMod10)

    }
  }


  complete_peak_pattern <- cp
  #complete_peak_pattern <- complete_peak_pattern[-which(round(complete_peak_pattern)==0)]
  #assign("complete_peak_pattern", complete_peak_pattern, envir = .GlobalEnv)

  span <- percentage_second_curve_height / 100
  if(span > 0.06) {
    span <- 0.1
  }
  df <- data.frame(x=1:length(complete_peak_pattern), y=complete_peak_pattern)
  loessMod <- loess(y ~ x, data=df, span=0.06, surface="direct") # 6% smoothing span
  # get smoothed output
  overlay <- predict(loessMod)

  if(length(which(overlay < 0)) > 0) {


    df <- data.frame(x=1:length(complete_peak_pattern), y=complete_peak_pattern)
    loessMod <- loess(y ~ x, data=df, span=0.04, surface="direct") # 6% smoothing span
    #overlay <- predict(loessMod)

    #overlay[which(overlay < 0)] <- 0

    df <- data.frame(x=1:length(overlay), y=overlay)
    loessMod <- loess(y ~ x, data=df, span=0.05, surface="direct") # 6% smoothing span
    # get smoothed output
    #overlay <- predict(loessMod)

  }


  #overlay <- complete_peak_pattern

  #assign("overlay", overlay, envir = .GlobalEnv)
  #assign("x_overlay_positions", x_overlay_positions, envir = .GlobalEnv)


  #y_second_curve <- y_second_curve * 0

  #}

  #goodness of fit
  if(alt.scores==T) {

    RSS <- sum(abs(heights_single_peaks - peak_heights_dtw_peaks) ^ 2)
    max_RSS <- round(sum((heights_single_peaks - c(0,0,0,0,0,0,0))^2))
    RSS_score <- 100 - ((RSS / max_RSS) * 100)

    #Perturbance value
    sample_proportion <- heights_sample_peaks / sum(heights_sample_peaks) * 100
    fit <- (dnorm(1:no.peaks, mean = 4, sd = 1.1)*10)*10
    perturbation_index <- sum(abs(fit-sample_proportion))
    max_perturbation <- sum(abs(fit-c(100,0,0,00,0,0,0)))
    normalized_perturbation_index <- 100 - ((perturbation_index / max_perturbation) * 100)

    #png(paste("~/Dropbox/LUmc/Genescan/WT_database/Selected_for_paper/Perturbation_VB17_3",".png", sep=""), width=150, height=150, units='mm', res=250)
    #plot(sample_proportion, xlim=c(1,7.3), ylim=c(0,40),type="h",lwd=16, xlab="peak number", ylab="proportion (%)", main="Perturbance Value")
    #points(c(1.3, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3),fit,type="h",lwd=16, col="red")

    #legend("topright",
    #       c(paste("Perturbance value :", round(perturbation_index,2)),
    #         paste("Normalized perturbance value :", round(normalized_perturbation_index,2)))
    #)

    #dev.off()

  }
  else {
    RSS_score <- 0
    perturbation_index <- 0
  }

  #######
  #calculate percentage difference between the peaks
  #peak_area_sample[which(peak_area_sample == 0 )] <- peak_area_sample[which(peak_area_sample == 0 )] + 0.1

  correction <- score_perfect_plot(three.bp.positions=36, no.peaks=no.peaks, percentage_second_curve_height=10, plot=F)

  peak_area_sample[which(peak_area_sample < (max(peak_area_sample) * 0.01))] <- max(peak_area_sample) * 0.01

  fraction <- peak_area_wt / (peak_area_sample)
  #wt peak area smaller then measured data area
  fraction[which(fraction < 1)] <- peak_area_wt[which(fraction < 1)] / peak_area_sample[which(fraction < 1)]
  #less penalty for peaks bigger then curve because this adds complexity
  fraction[which(fraction < 1)] <- (peak_area_sample[which(fraction < 1)] / peak_area_wt[which(fraction < 1)])

  if(is.infinite(max(fraction))) {
    #fraction[which(is.infinite(fraction))] <- sds[is.infinite(fraction)]
    fraction[which(is.infinite(fraction))] <- 100

  }

  #Fraction of difference of peak heights vs bell curve
  #peak_heights_dtw_peaks[which(peak_heights_dtw_peaks == 0 )] <- peak_heights_dtw_peaks[which(peak_heights_dtw_peaks == 0 )] + 0.1
  peak_heights_dtw_peaks[which(peak_heights_dtw_peaks < (max(peak_heights_dtw_peaks) * 0.01))] <- max(peak_heights_dtw_peaks) * 0.01

  fraction_heights <- peak_heights_dtw_peaks / heights_single_peaks
  #wt peak area smaller then measured data area
  fraction_heights[which(fraction_heights < 1)] <- peak_heights_dtw_peaks[which(fraction_heights < 1)] / heights_single_peaks[which(fraction_heights < 1)]
  #less penalty for peaks bigger then curve because this adds complexity
  fraction_heights[which(fraction_heights < 1)] <- (heights_single_peaks[which(fraction_heights < 1)] / peak_heights_dtw_peaks[which(fraction_heights < 1)])

  if(is.infinite(max(fraction_heights))) {
    #fraction_heights[which(is.infinite(fraction_heights))] <- sds[is.infinite(fraction_heights)]
    fraction_heights[which(is.infinite(fraction_heights))] <- 100
  }

  #fraction <- fraction * fraction_heights


  plot_area <-trapz(all_basepair_positions[[1]]$yy[pos_template][x])
  bell_curve_area <- trapz(y)
  second_curve_area <- trapz(y_second_curve)
  #second_curve_area <- 0

  #0.9360127
  #within_correction <- (plot_area - second_curve_area) / sum(peak_area_sample)
  within_correction <- correction[1]

  area_in_between_peaks<- (plot_area - second_curve_area) - (sum(peak_area_sample) * within_correction)
  #area_in_between_peaks <- area_in_between_peaks + sum(single_peak_penalties)
  area_in_between_wt_peaks <- (bell_curve_area - second_curve_area) - sum(peak_area_wt)

  percentage_in_between <- area_in_between_peaks / (area_in_between_wt_peaks / 100)


  plot_score <- (length(fraction) / sum(fraction)) * 100
  plot_score_heights <- (length(fraction_heights) / sum(fraction_heights)) * 100
  plot_score <- (plot_score + plot_score_heights)/2

  plot_score_correction <- correction[2]
  plot_score <- plot_score * plot_score_correction

  plot_score_in_between_penalty <- plot_score - (plot_score*((percentage_in_between*2)/100))


  if(plot==T) {
    points(x_overlay_positions[1]:(x_overlay_positions[1]+(length(overlay)-1)),
           overlay, type="l", col="firebrick", lwd=2)

    if (alt.scores == T) {
      legend("topright",
             c(paste("Peak score :", round(plot_score_in_between_penalty,2)),
               paste("GOF :", round(RSS_score,2)),
               paste("PV :", round(normalized_perturbation_index,2))
             ))

    }
    else {

      legend("topright",
             c(paste("Peak score :", round(plot_score_in_between_penalty,2))
             )
      )
    }
  }

  #peak_diff <- (heights_sample_peaks / heights_single_peaks)
  #peak_diff[which(peak_diff < 1)] <- (heights_single_peaks[which(peak_diff < 1)] / heights_sample_peaks[which(peak_diff < 1)])
  #peak_diff <- as.list((peak_diff))

  #names(peak_diff) <- c("1","2","3","4","5","6","7")
  #score_df <- data.frame(ps=plot_score_in_between_penalty , gof=RSS_score , pi=normalized_perturbation_index,peak_diff=peak_diff)


  #score_df <- data.frame(ps=plot_score_in_between_penalty)
  #score_df


  #plot(complete_peak_pattern,type="l")
  plot_score_in_between_penalty
}


