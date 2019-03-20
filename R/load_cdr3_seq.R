#' @title Load CDR3 sequencing data
#'
#' @description Load CDR3 sequencing data and prepare datatypes for scoring
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples load_cdr3_seq(folder, cdr3Length_column, geneFamily_column, sep="")
#'
#' @export load_cdr3_seq
load_cdr3_seq <- function(folder, cdr3Length_column, geneFamily_column, sep=sep) {

  #folder <- "~/Desktop/LUMC/CDR3_files/C57BL6_WT_set/"
  files <- list.files(folder)

  all_files <- c()

  pb <- txtProgressBar(min = 0, max = length(files), style = 3)

  for (i in 1:length(files)) {

    message("\nLoading file: ", files[i])

    if (i == 1) {

      cdr3 <- read.delim(paste(folder, files[i], sep = sep))
      cdr3 <- cdr3[,c(cdr3Length_column,geneFamily_column)]
      cdr3$sample <- gsub("\\..*","",files[i])

      all_files <- cdr3

    } else {

      cdr3 <- read.delim(paste(folder, files[i], sep = sep))
      cdr3 <- cdr3[,c(cdr3Length_column,geneFamily_column)]

      cdr3$sample <- gsub("\\..*","",files[i])

      all_files <- rbind(all_files,cdr3)

    }

    setTxtProgressBar(pb, i)

  }


  colnames(all_files)[3] <- "sample"

  families_seq <- all_files[,2]
  families_seq <- as.character(families_seq)
  all_files$families <- families_seq
  families_seq <- unique(families_seq)
  families_seq <- families_seq[order(families_seq)]

  samples <- unique(all_files$sample)

  seq_samples <- c()


  pb <- txtProgressBar(min = 0, max = length(families_seq), style = 3)
  message("\nCreating final data object")

  for (i in 1:length(families_seq)) {

    for (j in 1:length(samples)) {

      cdr3 <- all_files[which(all_files$sample == samples[j]),]
      table <- table(cdr3$cdr3Length[which(cdr3$families == families_seq[i])])

      if (length(table) > 1) {

        cdr3_df <- data.frame(cdr3_lengths=seq(as.numeric(names(table))[1], as.numeric(names(table))[length(table)], by=0.1), counts=NA)
        cdr3_df$counts[which(cdr3_df$cdr3_lengths %in% as.numeric(names(table)))] <- table

        interpolate <- approx(cdr3_df$counts, xout=seq_along(cdr3_df$counts))
        cdr3_df$counts[interpolate$x] <- interpolate$y


        if (max(table) > 400) {
          #LOESS
          loessMod <- loess(counts ~ cdr3_lengths, data=cdr3_df, span=0.03, surface="direct") # 6% smoothing span
          #get smoothed output
          pattern <- predict(loessMod)
          cdr3_df$counts <- pattern

          cdr3_df_big <- data.frame(cdr3_lengths=seq(0, 300, by=0.1), counts=NA)
          cdr3_df_big[which(as.character(cdr3_df_big$cdr3_lengths) %in% as.character(cdr3_df$cdr3_lengths)),] <- cdr3_df
          cdr3_df_big$counts[which(is.na(cdr3_df_big$count))] <- 0
          cdr3_df <- cdr3_df_big

          cdr3_list <- list(list(xx=cdr3_df$cdr3_lengths,yy=cdr3_df$counts))

          names(cdr3_list) <- paste(samples[j],"_", families_seq[i], sep="")
          seq_samples <- append(seq_samples,cdr3_list)

        } else{

          cdr3_df_big <- data.frame(cdr3_lengths=seq(0, 300, by=0.1), counts=NA)
          cdr3_df_big[which(as.character(cdr3_df_big$cdr3_lengths) %in% as.character(cdr3_df$cdr3_lengths)),] <- cdr3_df
          cdr3_df_big$counts[which(is.na(cdr3_df_big$count))] <- 0
          cdr3_df <- cdr3_df_big

          cdr3_list <- list(list(xx=cdr3_df$cdr3_lengths,yy=cdr3_df$counts))

          names(cdr3_list) <- paste(samples[j],"_", families_seq[i], sep="")
          seq_samples <- append(seq_samples,cdr3_list)
        }

      }

    }

    setTxtProgressBar(pb, i)
  }

  seq_samples

}
