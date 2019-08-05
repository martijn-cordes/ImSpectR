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
score_pattern_alignment <- function(bp_positions_dtw_peaks, no.peaks, three.bp.positions, total.pattern.width, plot.curve.fitting=T, plot.expected.model =T , plot=T, alt.scores=F)
{
  sd_correction <- c(1,1)
  if(length(all_basepair_positions[[1]]$xx) <= 3001) {
    sd_correction <- c(1.5,20)
  }

  sample <- all_basepair_positions[[1]]$yy[pos_template]

  score_curve_start <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)[1]
  score_curve_start <- score_curve_start - round(three.bp.positions / 1)
  if(score_curve_start < 0) {
    score_curve_start <- 1
  }
  score_curve_end <- (score_curve_start+total.pattern.width) + round(three.bp.positions / 1)

  if (all_basepair_positions[[1]]$xx[pos_template][score_curve_end] < bp_positions_dtw_peaks[length(bp_positions_dtw_peaks)]) {
    score_curve_end <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks[length(bp_positions_dtw_peaks)])+ round(three.bp.positions / 2)
  }

  if (score_curve_start-100 > 0) {
    x <- seq(score_curve_start-100,score_curve_end, by = 1)
  } else {
    x <- seq(score_curve_start-score_curve_start,score_curve_end, by = 1)
  }


  mean <- mean(x)
  mean <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)[4]

  sd <- total.pattern.width/5.3

  y <- dnorm(x, mean = mean, sd = sd)

  if( !is.na(max(sample) / max(sample[(score_curve_start+5):score_curve_end]) &&  max(sample) / max(sample[(score_curve_start+5):score_curve_end]) > 4)){
    y <- y * (max(sample[(score_curve_start+5):score_curve_end]) / max(y))
  }else {
    y <- y * (max(sample) / max(y))
  }

  if(plot.curve.fitting==T) {  points(x,y, type="l") }
  #if(plot==T) {  points(x,y, type="l") }

  position_peak <- which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks)
  max_peak_no <- which(y[x %in% position_peak] == max(y[x %in% position_peak]))

  peak_heights_dtw_peaks <- c()
  for (p in 1:length(bp_positions_dtw_peaks)) {
    peak_heights_dtw_peaks <- c(peak_heights_dtw_peaks,
                                all_basepair_positions[[1]]$yy[pos_template][which(all_basepair_positions[[1]]$xx[pos_template] %in% bp_positions_dtw_peaks[p] )]
    )
  }

  #Fit  all normal curves
  min_heights <- c()
  min_positions <- c()
  for(pp in 1:length(position_peak)) {
    min_position <- which(x %in% (position_peak[pp]-(three.bp.positions-10))) : which(x %in% (position_peak[pp]+(three.bp.positions+10)))
    min_positions <- c(min_positions, list(min_position))

    minimum <-min(all_basepair_positions[[1]]$yy[pos_template][(position_peak[pp]-(three.bp.positions-10)) : (position_peak[pp]+(three.bp.positions+10))])
    min_heights <- c(min_heights, minimum)
  }

  second_curve_height <- min_heights[round(length(min_heights)/2)] / max(peak_heights_dtw_peaks)

  y_second_curve <- dnorm(x, mean = mean, sd = total.pattern.width/4)
  y_second_curve <- y_second_curve * (((max(y) * second_curve_height) / max(y_second_curve)) *1.2)


  if((min_heights[2] < min(y_second_curve[min_positions[[2]] ]))) {
    y_second_curve <- y_second_curve / (min(y_second_curve[min_positions[[2]]]) / min_heights[2])
  }

  if ((max(y_second_curve) / max(y)) > 0.3 | (min_heights[2] < min(y_second_curve[min_positions[[2]] ]))) {
    y_second_curve <- y_second_curve * (((max(y) * 0.05) / max(y_second_curve)) *1.5)
  }

  y_second_curve <- y_second_curve-min(y_second_curve)

  if(plot.curve.fitting == T) {
    points(x,y_second_curve, type="l", col="red")
  }


  #Start scoring

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

    x_single_peak <- seq(mean_single_peak-20,mean_single_peak+20, by = 1)

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
    if(percentage_second_curve_height < 1) {
      percentage_second_curve_height <- 5
    }

    sd <- round(three.bp.positions / (10-(log2(ceil(percentage_second_curve_height) +1)) ))
    slope_heightest_peak <- (y[mean(which(round(x) == position_peak[max_peak_no][1])[1])] / sd)
    sd <- y[mean(which(round(x) == position_peak[p])[1])] / slope_heightest_peak
    sd <- sd*sd_correction[1] + (log2((max(y)/ y[mean(which(round(x) == position_peak[p])[1])]))*(log10(ceil(percentage_second_curve_height) + sd_correction[2])))
    sds <- c(sds, sd)

    height_single_peak <- y[mean(which(round(x) == position_peak[p])[1])]

    y_single_peak <- dnorm(x_single_peak, mean = mean_single_peak, sd = sd)
    y_single_peak  <- y_single_peak  * height_single_peak / max(y_single_peak )

    overlaps <- order(abs(y_second_curve[which(x %in% x_single_peak)] - y_single_peak))
    overlaps <- overlaps[c(1,(which(abs(diff(overlaps)) > 5) + 1)[1])]

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

    if(plot==T){
    # abline(v=x_single_peak[1],lty="dotted")
    # abline(v=x_single_peak[length(x_single_peak)], lty="dotted")


    # polygon(x = c(x_single_peak[1],x_single_peak, x_single_peak[length(x_single_peak)]),
    #        y = c(y_second_curve[which(x %in% x_single_peak[1])],
    #            y_single_peak,
    #            y_second_curve[which(x %in% x_single_peak[length(x_single_peak)])] ),
    #        col = "deepskyblue2", border="deepskyblue2")
    #
    #
    # polygon(x = c(x_single_peak[1], x_single_peak[1:(length(x_single_peak))], x_single_peak[length(x_single_peak)]),
    #        y = c(y_second_curve[which(x %in% x_single_peak[1])],
    #              all_basepair_positions[[1]]$yy[pos_template][x_single_peak[1]:x_single_peak[(length(x_single_peak))]],
    #              y_second_curve[which(x %in% x_single_peak[length(x_single_peak)])] ),
    #      col = "chartreuse2", border="chartreuse2")



    }

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
    peak_area_wt <- c(peak_area_wt, trapz(y_single_peak))
    peak_area_sample <- c(peak_area_sample, peak_area_corrected_sample)
    peak_area_second_curve_single_peaks <- c(peak_area_second_curve_single_peaks, trapz(y_second_curve[which(x %in% x_single_peak)]))

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

  for (g in 1:length(groups)) {

    if(g == 1) {
      data <- cp[(unlist(groups[g])[1]) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct")
      cp[(unlist(groups[g])[1]) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)] <- predict(loessMod10)
    }

    else if(g == length(groups)) {
      data <- cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))])]
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct")
      cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))])] <- predict(loessMod10)
    }

    else {
      data <- cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]
      data <- data.frame(x=1:length(data), y=data)
      loessMod10 <- loess(y ~ x, data=data, span=1, surface="direct")
      cp[(unlist(groups[g])[1]-5) : (unlist(groups[g])[length(unlist(groups[g]))] + 5)]  <- predict(loessMod10)
    }
  }


  complete_peak_pattern <- cp

  span <- percentage_second_curve_height / 100
  if(span > 0.06) {
    span <- 0.1
  }
  df <- data.frame(x=1:length(complete_peak_pattern), y=complete_peak_pattern)
  loessMod <- loess(y ~ x, data=df, span=0.06, surface="direct") # 6% smoothing span
  # get smoothed output
  overlay <- predict(loessMod)

  if(length(which(overlay < 0)) > 0) {
    overlay[which(overlay < 0)] <- 0
  }


  #goodness of fit
  if(alt.scores==T) {

    RSS <- sum(abs(heights_single_peaks - peak_heights_dtw_peaks) ^ 2)
    max_RSS <- round(sum((heights_single_peaks - c(0,0,0,0,0,0,0))^2))
    RSS_score <- 100 - ((RSS / max_RSS) * 100)

    chisq <- sum((as.numeric(peak_heights_dtw_peaks)-heights_single_peaks)^2/heights_single_peaks)
    max_chisq <- sum( (as.numeric(0)-heights_single_peaks)^2/heights_single_peaks)
    chisq_score <- 100 - (chisq/(max_chisq) * 100)


    #(value-min)/(max-min)

    #Perturbance value
    sample_proportion <- heights_sample_peaks / sum(heights_sample_peaks) * 100
    fit <- (dnorm(1:no.peaks, mean = 4, sd = 1.1)*10)*10
    perturbation_index <- sum(abs(fit-sample_proportion))
    max_perturbation <- sum(abs(fit-c(100,0,0,00,0,0,0)))
    normalized_perturbation_index <- 100 - ((perturbation_index / max_perturbation) * 100)

  }
  else {
    chisq_score <- 0
    perturbation_index <- 0
  }

  #######
  #calculate percentage difference between the peaks
  correction <- score_perfect_plot(three.bp.positions=36, no.peaks=no.peaks, percentage_second_curve_height=10, plot=F)

  peak_area_sample[which(peak_area_sample < (max(peak_area_sample) * 0.01))] <- max(peak_area_sample) * 0.01

  fraction <- peak_area_wt / (peak_area_sample)
  #wt peak area smaller then measured data area
  fraction[which(fraction < 1)] <- peak_area_wt[which(fraction < 1)] / peak_area_sample[which(fraction < 1)]
  #less penalty for peaks bigger then curve because this adds complexity
  fraction[which(fraction < 1)] <- (peak_area_sample[which(fraction < 1)] / peak_area_wt[which(fraction < 1)])

  if(is.infinite(max(fraction))) {
    fraction[which(is.infinite(fraction))] <- 100
  }

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
  all_basepair_positions[[1]]$yy[pos_template][x] <- all_basepair_positions[[1]]$yy[pos_template][x] - y_second_curve
  all_basepair_positions[[1]]$yy[pos_template][x][which(all_basepair_positions[[1]]$yy[pos_template][x] < 0)] <- 0
  y <- y - y_second_curve


  ####Calculate areas between each peak and maximum area possible between each peak.
  in_between_seperate <- c()
  in_between_max <- c()

  in_between_seperate <- c(in_between_seperate, sum(all_basepair_positions[[1]]$yy[pos_template][(all_overlaps[[1]][1]-20):all_overlaps[[1]][1]]))
  in_between_max <- c(in_between_max, sum(y[(which(x == all_overlaps[[1]][1])-20) : which(x == all_overlaps[[1]][1])]))


  #####PLOT ILLUSTRATIONS
  original_plot <- all_basepair_positions[[1]]
  original_plot$yy[pos_template][x] <- original_plot$yy[pos_template][x] + y_second_curve
  y_original <- y + y_second_curve


  if(plot==T) {

    # polygon(x = c((all_overlaps[[1]][1]-20),(all_overlaps[[1]][1]-20):(all_overlaps[[1]][1]),all_overlaps[[1]][1] ),
    #         y = c(y_second_curve[which(x %in% (all_overlaps[[1]][1]-20))],
    #               y_original[which(x %in% (all_overlaps[[1]][1]-20):(all_overlaps[[1]][1]-0))],
    #               y_second_curve[which(x %in% all_overlaps[[1]][1])]),
    #         col = "royalblue4", border="royalblue4" )
    #
    # polygon(x = c((all_overlaps[[1]][1]-20):(all_overlaps[[1]][1]),all_overlaps[[1]][1] ),
    #         y = c(original_plot$yy[pos_template][(all_overlaps[[1]][1]-20):(all_overlaps[[1]][1]-0)],
    #               y_second_curve[which(x %in% all_overlaps[[1]][1])]),
    #         col = "firebrick1", border="firebrick1")
    #
  }
  ###########

  for(o in 1:(length(all_overlaps)-1)) {
    in_between_seperate <- c(in_between_seperate, sum(all_basepair_positions[[1]]$yy[pos_template][all_overlaps[[o]][2]:all_overlaps[[o+1]][1]]))
    in_between_max <- c(in_between_max, sum(y[which(x == all_overlaps[[o]][2]) : which(x == all_overlaps[[o+1]][1])]))

    #################
    if (plot==T) {

      # polygon(x = c(all_overlaps[[o]][2],(all_overlaps[[o]][2]:all_overlaps[[o+1]][1]),all_overlaps[[o+1]][1]),
      #        y = c(y_second_curve[which(x %in% all_overlaps[[o]][2])],
      #              y_original[which(x %in% (all_overlaps[[o]][2]):all_overlaps[[o+1]][1])],
      #              y_second_curve[which(x %in% all_overlaps[[o+1]][1])]),
      #        col = "royalblue4", border="royalblue4" )
      #
      #
      # polygon(x = c(all_overlaps[[o]][2],(all_overlaps[[o]][2]:all_overlaps[[o+1]][1]),all_overlaps[[o+1]][1]),
      #         y = c(y_second_curve[which(x %in% all_overlaps[[o]][2])],
      #               original_plot$yy[pos_template][(all_overlaps[[o]][2]):all_overlaps[[o+1]][1]],
      #               y_second_curve[which(x %in% all_overlaps[[o+1]][1])]),
      #         col = "firebrick1", border="firebrick1")


    }
    ##########################


}

  in_between_seperate <- c(in_between_seperate, sum(all_basepair_positions[[1]]$yy[pos_template][all_overlaps[[length(all_overlaps)]][2]:(all_overlaps[[length(all_overlaps)]][2]+20)]))
  in_between_max <- c(in_between_max, sum(y[which(x == all_overlaps[[length(all_overlaps)]][2]) : (which(x == all_overlaps[[length(all_overlaps)]][2]+20))]))

  #######
  if (plot==T)   {
    # polygon(x = c(all_overlaps[[length(all_overlaps)]][2],(all_overlaps[[length(all_overlaps)]][2]:(all_overlaps[[length(all_overlaps)]][2]+20)),(all_overlaps[[length(all_overlaps)]][2]+20)),
    #         y = c(y_second_curve[which(x %in% all_overlaps[[length(all_overlaps)]][2])   ],
    #               y_original[which(x %in% (all_overlaps[[length(all_overlaps)]][2]:(all_overlaps[[length(all_overlaps)]][2]+20)))],
    #               y_second_curve[which(x %in% (all_overlaps[[length(all_overlaps)]][2]+20))]),
    #         col = "royalblue4", border="royalblue4" )
    #
    # polygon(x = c(all_overlaps[[length(all_overlaps)]][2],(all_overlaps[[length(all_overlaps)]][2]:(all_overlaps[[length(all_overlaps)]][2]+20)),(all_overlaps[[length(all_overlaps)]][2]+20)),
    #         y = c(y_second_curve[which(x %in% all_overlaps[[length(all_overlaps)]][2])   ],
    #               original_plot$yy[pos_template][(all_overlaps[[length(all_overlaps)]][2]:(all_overlaps[[length(all_overlaps)]][2]+20))],
    #               y_second_curve[which(x %in% (all_overlaps[[length(all_overlaps)]][2]+20))]),
    #         col = "firebrick1", border="firebrick1")

  }

  ###########
  sample_data <- all_basepair_positions[[1]]$yy[pos_template][x_overlay_positions[1]:(x_overlay_positions[1]+(length(overlay)-1))]

  ###Score calculation
  plot_score <- (length(fraction) / sum(fraction)) * 100
  score_1 <- plot_score
  score_1 <- plot_score * correction[1]

  plot_score_heights <- (length(fraction_heights) / sum(fraction_heights)) * 100
  plot_score_heights <- plot_score_heights * correction[2]


  residual_penalty <- (trapz(in_between_seperate) / trapz(in_between_max)*100)
  residual_penalty <- abs(residual_penalty-correction[3])

  combined_score <- ((score_1+plot_score_heights)/2) - residual_penalty
  combined_score <- score_1 - residual_penalty

  if(plot.expected.model==T) {
    points(x_overlay_positions[1]:(x_overlay_positions[1]+(length(overlay)-1)),
           overlay, type="l", col="firebrick", lwd=2, lty="dashed")
  }


  if(plot==T) {
    if (alt.scores == T) {
      # legend("topright",
      #        c(paste("ImSpectR score :", round(combined_score,2)),
      #          paste("GOF :", round(chisq_score,2)),
      #          paste("PV :", round(normalized_perturbation_index,2))
      #        ))

    }
    else {
       legend("topleft",legend=c("Sample", "Expected Model"),
             col=c("black", "red"),lty=1, cex=0.8)

       legend("topright",
              c( paste("Peak score:", round(combined_score ,2))
              ))


      # legend("topright",
      #       c(paste("Global pattern score:", round(plot_score_heights,2)),
      #         paste("Individual peak score:", round(score_1,2)),
      #         paste("Residual penalty:", round(residual_penalty,2)),
      #         paste("Peak score:", round(combined_score ,2)),
      #         paste("Peak score without overall:", round(score_1-residual_penalty ,2))
      #
      #       )
      # )
    }
  }

  #if (plot == T) {
  #  points(sample , col="black" ,
  #         type="l",
  #         xlim=c(x_positions_dtw_peaks[1]-100, x_positions_dtw_peaks[length(x_positions_dtw_peaks)]+100), lwd=2)
  #}


  if(alt.scores == T) {
    list(peak_score = combined_score,GOF = chisq_score,PI=normalized_perturbation_index)

  } else {
    combined_score
  }

}


