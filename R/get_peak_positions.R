#' @title Get peak positions
#'
#' @description Get peak positions of the peaks matched by Dynamic Time Warping
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples get_peak_positions(alignment, sample)
#'
#' @export get_peak_positions
get_peak_positions <- function(alignment, sample, peak.margin=NULL,plot=T,plot.curve.fitting=plot.curve.fitting, plot.expected.model=T)
{
  if (plot==F) {
    plot.expected.model <- F
  }


  if (is.null(peak.margin)) {
    peak.margin <- 10
  }

  all_basepair_positions <- sample
  assign("all_basepair_positions", all_basepair_positions, envir = .GlobalEnv)


  if(is.null(pos_template)){
    assign("pos_template", (1:length(all_basepair_positions[[1]]$xx)), envir = .GlobalEnv)
  }


  sample <- all_basepair_positions[[1]]$yy[alignment$x_pos]

  #x coordinates for peaks of dtw matched pattern
  peaks <- fpeaks(alignment$y, f=200, amp=c(100,100), plot=F)

  if (nrow(peaks) < no.peaks){
    peaks <- fpeaks(alignment$y, f=200, amp=c(40,0.01), plot=F)
  }

  if (nrow(peaks) < no.peaks){
    peaks <- fpeaks(alignment$y, f=200, amp=c(10,10), plot=F)
  }

  if (nrow(peaks) < no.peaks){
    peaks <- fpeaks(alignment$y, f=200, amp=c(1,1), plot=F)
  }

  if (nrow(peaks) < no.peaks){
    peaks <- fpeaks(alignment$y, f=200, amp=c(0.1,0.1), plot=F)
  }

  if (nrow(peaks) < no.peaks){
    peaks <- fpeaks(alignment$y, f=200, amp=c(0.01,0.01), plot=F)
  }

  bp_positions_dtw_peaks <- c()
  x_positions_dtw_peaks <- c()
  peak_heights_dtw_peaks <- c()
  duplicated <- which(duplicated(peaks[,2]))
  if (isempty(duplicated)) {
    duplicated <- 0
  }

  for (p in 1:nrow(peaks)) {
    if (p == duplicated) {
      position <- alignment$x[which(alignment$y %in% peaks[,2][p])][2]
    } else {
      position <- alignment$x[which(alignment$y %in% peaks[,2][p])][1]
    }

    if (is.na(position))  {
      position <- alignment$x[which(floor(alignment$y) %in% floor(peaks[,2][p]))][1]
    }
    if (is.na(position))  {
      position <- alignment$x[which(floor(alignment$y) %in% round(peaks[,2][p],-1))][1]
    }

    if (position - peak.margin <= 0) {
      peak.margin <- position
    }

    if (position + peak.margin > length(sample)) {
      peak.margin <- length(sample) - position
    }


    if (max(sample[(position-peak.margin) : (position+peak.margin)], na.rm=T) > 0 ) {
      x_pos <- which(all_basepair_positions[[1]]$yy[pos_template] %in% sample[which(sample == max(sample[(position-peak.margin) : (position+peak.margin)], na.rm=T))])

      x_positions_dtw_peaks <- c(x_positions_dtw_peaks,x_pos)

      peak_heights_dtw_peaks <- c(peak_heights_dtw_peaks,
                                  sample[which(sample == max(sample[(position-peak.margin) : (position+peak.margin)], na.rm=T))]
      )

      bp_positions_dtw_peaks <- c(bp_positions_dtw_peaks,
                                  all_basepair_positions[[1]]$xx[pos_template][which(all_basepair_positions[[1]]$yy[pos_template] %in% sample[which(sample == max(sample[(position-peak.margin) : (position+peak.margin)], na.rm=T))])]
      )

    }
    else {
      x_pos <- which(pos_template %in% alignment$x_pos[which(alignment$y %in% peaks[,2][p])])
      x_positions_dtw_peaks <- c(x_positions_dtw_peaks,x_pos)

      peak_heights_dtw_peaks <- c(peak_heights_dtw_peaks, 0)

      bp_positions_dtw_peaks <- c(bp_positions_dtw_peaks,
                                  alignment$bp[which(alignment$y %in% peaks[,2][p])])
    }
  }


  if(plot == T) {

    plot(all_basepair_positions[[1]]$yy[pos_template] , col="black" , type="l", main=names(all_basepair_positions)[1], xlim=c(x_positions_dtw_peaks[1]-100, x_positions_dtw_peaks[length(x_positions_dtw_peaks)]+100),ylim=c(0,max(sample)), xaxt="n", xlab="length in basepairs", ylab="", lwd=2)
    text(x_positions_dtw_peaks, par("usr")[3] - 0.2, labels =round(bp_positions_dtw_peaks) , srt = 45, pos = 1, xpd = TRUE)

    if (plot.curve.fitting == T) {
      points(
        x_positions_dtw_peaks,
        peak_heights_dtw_peaks, col="red"
      )
    }

  }

  bp_positions_dtw_peaks
}
