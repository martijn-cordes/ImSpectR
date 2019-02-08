#' @title Score a complete dataset
#'
#' @description Score a complete dataset by iterating score_sample. See score_sample for functional output.
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples score_dataset(dataset, no.peaks)
#'
#' @export score_dataset
score_dataset <- function(dataset, no.peaks, peak.margin=NULL, peak.window=NULL, window.size=NULL, plot.pattern.matching =F, plot.curve.fitting=F, plot.expected.model=T, plot=T ) {


  scored_dataset <- c()
  pb <- txtProgressBar(min = 0, max = length(dataset), style = 3)
  for (i in 1:length(dataset)) {
    message("\nAnalyzing sample: ",names(dataset[i]))


    scored_sample <- score_sample(dataset[i], no.peaks, peak.margin=peak.margin, peak.window=peak.window, window.size=window.size, plot.pattern.matching = plot.pattern.matching, plot.curve.fitting = plot.curve.fitting, plot.expected.model = plot.expected.model, plot=plot)
    scored_dataset <- rbind(scored_dataset, scored_sample)
    setTxtProgressBar(pb, i)
  }

  scored_dataset
}

