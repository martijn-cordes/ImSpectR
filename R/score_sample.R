#' @title Score a single sample
#'
#' @description Score a single sample. This function will search the data for a peak pattern best fitting an in-silico pattern defined by the user, select the in-frame peaks and score the pattern.
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples score_sample(sample, no.peaks)
#'
#' @export score_sample
score_sample <- function(sample, no.peaks, peak.margin=NULL, peak.window=NULL, window.size=NULL, plot.pattern.matching=F, plot.curve.fitting=F, plot=T) {

  if(length(sample)>1) {
    current_sample <- score_dataset(sample, no.peaks, peak.margin=NULL, peak.window=NULL, window.size=NULL, plot.pattern.matching = F, plot.curve.fitting = plot.curve.fitting, plot=plot)
  }

  else {

    if(is.null(peak.margin)) {
      peak.margin <- 10
    }

    if(is.null(peak.window)) {
      peak.window <- c(100,300)
    }

    if(is.null(window.size)) {
      window.size <- 1
    }


    all_basepair_positions <- sample
    sample_name <- names(all_basepair_positions)

    three.bp.positions <- length(which((all_basepair_positions[[1]]$xx > peak.window[1]) & (all_basepair_positions[[1]]$xx < (peak.window[1]+3))))
    total.pattern.width <- three.bp.positions * no.peaks

    #Select window to search the pattern in from the whole dataset
    position_template <- which((all_basepair_positions[[1]]$xx > peak.window[1]) & (all_basepair_positions[[1]]$xx < peak.window[2]))
    assign("pos_template", position_template, envir = .GlobalEnv)


    sample_data <- all_basepair_positions[[1]]$yy[position_template]
    names(sample_data) <- sample_name

    if(max(sample_data) < 1000 | (three.bp.positions < 20) ) {
      message("No peak pattern found")

      peak_score <- 0
      #current_sample <- list(rawdata = all_basepair_positions, plot_data = all_basepair_positions, pattern_alignment=NA, peak_positions = NA, score=peak_score)
      #current_sample <- list(current_sample)
      #names(current_sample) <- sample_name
      current_sample <- data.frame(peak.score=peak_score)
      rownames(current_sample) <- sample_name

      if(plot==T) {
        plot( all_basepair_positions[[1]]$xx,  all_basepair_positions[[1]]$yy, type="l", xaxt="n",xlab="position in basepairs", ylab="" ,main=sample_name)
        axis(side=1,at=seq(0,length(all_basepair_positions[[1]]$xx),10))
      }

    } else {
      #Use a sliding window of the width of the pattern to match to slide through the dataset, performing DTW on every window
      #Two outcome values of DTW will be used to select pattern to use for scoring:
      #1-Width of matched pattern, the closer to the actual width of the pattern the better
      #2-Distance, which is the distance of the matrix value, the lowest distance is usually the best match.

      artifical_WT <- generate_artifical_WT(no.peaks, plot=F, three.bp.positions = three.bp.positions)
      artifical_WT <- artifical_WT * (max(all_basepair_positions[[1]]$yy[position_template] ) / max(artifical_WT))
      query <- artifical_WT

      #sliding window for matched pattern
      start <- peak.window[1]
      min.width.pattern <- (no.peaks * 4) + 2
      #min.width.pattern <- 30
      all_widths <- c()
      all_starts <- c()
      all_distances <- c()

      #total <- peak.window[2]-min.width.pattern
      #message("Performing sliding window DTW")
      #pb <- txtProgressBar(min = start, max = total, style = 3)

      while(start < peak.window[2]-min.width.pattern) {

        window_start <- which(all_basepair_positions[[1]]$xx > start)[1]
        window_end <- which(all_basepair_positions[[1]]$xx > start+min.width.pattern)[1]
        sample_data <- all_basepair_positions[[1]]$yy[window_start:window_end]

        if (max(sample_data) > 1000) {
          #pattern_alignment <- align_peak_pattern(query, sample_data, start=start, pattern.width = min.width.pattern, window.size = window.size, plot=F)
          pattern_alignment <- align_peak_pattern(query, all_basepair_positions, start=start, pattern.width = min.width.pattern, window.size = window.size,  plot=F)

          matched_pattern_range <- all_basepair_positions[[1]]$xx[pattern_alignment$x + which(all_basepair_positions[[1]]$xx > start)[1]]
          matched_pattern.width <- matched_pattern_range[length(matched_pattern_range)] - matched_pattern_range[1]

          all_widths <- c(all_widths,matched_pattern.width)
          all_distances <- c(all_distances, pattern_alignment$distance[1])

          start <- start + 1
          all_starts <- c(all_starts, start)



        } else {
          start <- start + 1
          all_starts <- c(all_starts, start)
          matched_pattern.width <- 0
          all_widths <- c(all_widths,matched_pattern.width)
          distance <- 100000
          all_distances <- c(all_distances, distance)
        }

      }


      #####
      #Get best match
      #####

      #A minimal width of a pattern has to be used
      if (max(all_widths)  < 10) {
        message("No peak pattern found")

        peak_score <- 0
        #current_sample <- list(rawdata = all_basepair_positions, plot_data = all_basepair_positions, pattern_alignment=NA, peak_positions = NA, score=peak_score)
        #current_sample <- list(current_sample)
        #names(current_sample) <- sample_name
        current_sample <- data.frame(peak.score=peak_score)
        rownames(current_sample) <- sample_name

        if(plot==T) {
          plot( all_basepair_positions[[1]]$xx,  all_basepair_positions[[1]]$yy, type="l", xaxt="n",xlab="position in basepairs", ylab="" ,main=sample_name)
          axis(side=1,at=seq(0,length(all_basepair_positions[[1]]$xx),10))
        }


      } else {
        start <- (peak.window[1] -1) + which(all_distances == min(all_distances))[1]

        #Which is the best mean between widths and distance
        #score_widths <- all_widths
        #score_widths[rev(order(score_widths))] <- 1 : length(all_widths)
        score_distances <- all_distances
        score_distances[order(score_distances)] <- 1 : length(all_widths)
        #scores <- (score_widths + score_distances)
        #start <- (peak.window[1] -1) + order(scores)[1]

        top_distances <- order(score_distances)[1:5]

        top_scores <- c()
        for (t in 1:length(top_distances)) {
          start <- (peak.window[1] -1) + top_distances[t]

          if(start <= (peak.window[1]+1)) {
            top_scores <- c(top_scores, 0)
          }

          else {

            window_start <- which(all_basepair_positions[[1]]$xx > start)[1]
            window_end <- which(all_basepair_positions[[1]]$xx > start+min.width.pattern)[1]
            sample_data <- all_basepair_positions[[1]]$yy[window_start:window_end]

            #pattern_alignment <- align_peak_pattern(query, sample_data, start=start, pattern.width=min.width.pattern, window.size= window.size, plot=F)
            pattern_alignment <- align_peak_pattern(query, all_basepair_positions, start=start, pattern.width = min.width.pattern, window.size = window.size,  plot=F)

            #bp_positions_dtw_peaks <- get_peak_positions(pattern_alignment, sample_data, peak.margin = peak.margin, plot=F)
            bp_positions_dtw_peaks <- get_peak_positions(pattern_alignment, all_basepair_positions, peak.margin = peak.margin, plot=F)


            if (nrow(pattern_alignment) > three.bp.positions * (no.peaks + 1)) {
              total.pattern.width <- nrow(pattern_alignment)
            } else {
              total.pattern.width <- three.bp.positions * (no.peaks + 1)
            }


            pm <- peak.margin
            while (min(diff(bp_positions_dtw_peaks)) < 1 | max(diff(bp_positions_dtw_peaks)) > 4) {
              pm <- pm - 1
              bp_positions_dtw_peaks <- get_peak_positions(pattern_alignment, all_basepair_positions, peak.margin = peak.margin, plot=F)

              if(pm ==0 ){
                break
              }

            }

            peak_score <- score_pattern_alignment(bp_positions_dtw_peaks, no.peaks, three.bp.positions, total.pattern.width, alt.scores = F, plot.curve.fitting = F, plot=F)
            top_scores <- c(top_scores, peak_score)
            # else {
            #  top_scores <- c(top_scores, 0)
            #}
          }

        }

        if (length(which(is.na(top_scores)) > 0)) {
          top_scores[which(is.na(top_scores))] <- 0
        }

        if(max(top_scores[1:3]) > 0) {
          top_scores <- top_scores[1:5]
        }


        #calculate final score of pattern with the best distance and the best score
        if (max(top_scores) > 0) {

          start <- (peak.window[1] -1) + top_distances[which(top_scores == max(top_scores))[1]]
          #start <- (peak.window[1] -1) + top_distances[4]

          #sample_starts <- c(sample_starts,start)

          window_start <- which(all_basepair_positions[[1]]$xx > start)[1]
          window_end <- which(all_basepair_positions[[1]]$xx > start+min.width.pattern)[1]
          sample_data <- all_basepair_positions[[1]]$yy[window_start:window_end]

          #png(paste("~/Dropbox/LUmc/Genescan/WT_database/Selected_for_paper/WTVB17_2",".png", sep=""), width=200, height=150, units='mm', res=250)
          pattern_alignment <- align_peak_pattern(query, all_basepair_positions, start=start, pattern.width=min.width.pattern, window.size= window.size, plot=plot.pattern.matching)
          bp_positions_dtw_peaks <- get_peak_positions(pattern_alignment, all_basepair_positions, peak.margin = peak.margin, plot=plot)
          #points(df_smooth$x + window_start,df_smooth$y,type="l", col="darkgreen", lwd=2)

          if (nrow(pattern_alignment) > three.bp.positions * (no.peaks + 1)) {
            total.pattern.width <- nrow(pattern_alignment)
          } else {
            total.pattern.width <- three.bp.positions * (no.peaks + 1)
          }

          peak_score <- score_pattern_alignment(bp_positions_dtw_peaks, no.peaks, three.bp.positions,total.pattern.width, plot.curve.fitting=plot.curve.fitting, plot=T, alt.scores=F)
          #dev.off()

          plot_score <- peak_score
          message("Scoring complete")
          #peak_score_df <- rbind(peak_score_df, c(fsa.names[i],as.numeric(peak_score)))

        } else{
          peak_score <- 0
          #current_sample <- list(rawdata = all_basepair_positions, plot_data = all_basepair_positions, pattern_alignment=NA, peak_positions = NA, score=peak_score)
          #current_sample <- list(current_sample)
          #names(current_sample) <- sample_name
          current_sample <- data.frame(peak.score=peak_score)
          #rownames(current_sample) <- sample_name

          if(plot==T) {
            plot( all_basepair_positions[[1]]$xx,  all_basepair_positions[[1]]$yy, type="l", xaxt="n",xlab="position in basepairs", ylab="" ,main=sample_name)
            axis(side=1,at=seq(0,length(all_basepair_positions[[1]]$xx),10))
          }

        }

        #current_sample <- list(rawdata = all_basepair_positions, plot_data = sample_data, pattern_alignment=pattern_alignment, peak_positions = bp_positions_dtw_peaks, score=peak_score$ps)
        #current_sample <- list(current_sample)
        #names(current_sample) <- sample_name
        current_sample <- data.frame(peak.score=peak_score)
        rownames(current_sample) <- sample_name

      }
    }

    pos_template <- NULL

    current_sample

  }

}
