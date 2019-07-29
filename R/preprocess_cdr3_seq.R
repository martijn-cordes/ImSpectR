#' @title Preprocess CDR3 sequencing data
#'
#' @description Preprocess CDR3 sequencing data and prepare datatypes for scoring
#'
#' @param symbol
#'
#' @return NULL
#'
#' @examples preprocess_cdr3_seq(folder, cdr3Length_column, geneFamily_column, sep="")
#'
#' @export preprocess_cdr3_seq
preprocess_cdr3_seq <- function(folder, cdr3Length_column, geneFamily_column, sep=NULL) {
  if(is.null(sep)) {
    sep <- "\t"
  }
  
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)

  if (substr(folder, nchar(folder),nchar(folder)) != "/") {
    folder <- paste(folder,"/",sep="")
  }
  
  files <- list.files(folder)

  all_files <- c()

  pb <- txtProgressBar(min = 0, max = length(files), style = 3)

  for (i in 1:length(files)) {

    message("\nLoading file: ", files[i])

    if (i == 1) {

      cdr3 <- read.delim(paste(folder, files[i], sep = ""), sep=sep)
      cdr3 <- cdr3[,c(cdr3Length_column,geneFamily_column)]
      cdr3$sample <- gsub("\\..*","",files[i])

      all_files <- cdr3

    } else {

      cdr3 <- read.delim(paste(folder, files[i], sep = ""), sep=sep)
      cdr3 <- cdr3[,c(cdr3Length_column,geneFamily_column)]

      cdr3$sample <- gsub("\\..*","",files[i])

      all_files <- rbind(all_files,cdr3)

    }

    setTxtProgressBar(pb, i)

  }


  colnames(all_files) <- c("cdr3Length","geneFamily", "sample")

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
        density_cdr3 <- density(cdr3$cdr3Length[which(cdr3$families == families_seq[i])],bw=0.5)
        
        cdr3_df <- data.frame(cdr3_lengths=seq(min(round(density_cdr3$x,2)), max(round(density_cdr3$x,2)), by=0.01), counts=NA)
        cdr3_df$counts[which(cdr3_df$cdr3_lengths %in% round(density_cdr3$x,2))] <- 
          density_cdr3$y[which(round(density_cdr3$x,2) %in% cdr3_df$cdr3_lengths)]
        
        #cdr3_df <- data.frame(cdr3_lengths=seq(as.numeric(names(table))[1], as.numeric(names(table))[length(table)], by=0.01), counts=NA)
        #cdr3_df$counts[which(cdr3_df$cdr3_lengths %in% as.numeric(names(table)))] <- table
        
        interpolate <- approx(cdr3_df$counts, xout=seq_along(cdr3_df$counts))
        cdr3_df$counts[interpolate$x] <- interpolate$y
        
        cdr3_df$counts <- cdr3_df$counts * max(table)/max(cdr3_df$counts, na.rm=T)
        
        cdr3_df <- cdr3_df[which(cdr3_df$cdr3_lengths %in% seq(min(round(density_cdr3$x,1)), max(round(density_cdr3$x,1)), by=0.1)),] 
        
        
        if (max(table) > 400) {
          #LOESS
          #loessMod <- loess(counts ~ cdr3_lengths, data=cdr3_df, span=0.03, surface="direct") # 6% smoothing span
          #get smoothed output
          #pattern <- predict(loessMod)
          #cdr3_df$counts <- pattern

          cdr3_df_big <- data.frame(cdr3_lengths=seq(0, 300, by=0.1), counts=NA)
          cdr3_df_big[which(as.character(cdr3_df_big$cdr3_lengths) %in% as.character(cdr3_df$cdr3_lengths)),] <- cdr3_df
          #cdr3_df_big$counts[which(is.na(cdr3_df_big$count))] <- 0
          cdr3_df <- cdr3_df_big
          
          interpolate <- approx(cdr3_df$counts, xout=seq_along(cdr3_df$counts))
          cdr3_df$counts[interpolate$x] <- interpolate$y
          cdr3_df$counts[which(is.na(cdr3_df$count))] <- 0
          
          cdr3_list <- list(list(xx=cdr3_df$cdr3_lengths,yy=cdr3_df$counts))
          
          names(cdr3_list) <- paste(samples[j],"_", families_seq[i], sep="")
          seq_samples <- append(seq_samples,cdr3_list)

        } else{

          cdr3_df_big <- data.frame(cdr3_lengths=seq(0, 300, by=0.1), counts=NA)
          cdr3_df_big[which(as.character(cdr3_df_big$cdr3_lengths) %in% as.character(cdr3_df$cdr3_lengths)),] <- cdr3_df
          #cdr3_df_big$counts[which(is.na(cdr3_df_big$count))] <- 0
          cdr3_df <- cdr3_df_big
          
          interpolate <- approx(cdr3_df$counts, xout=seq_along(cdr3_df$counts))
          cdr3_df$counts[interpolate$x] <- interpolate$y
          cdr3_df$counts[which(is.na(cdr3_df$count))] <- 0
          
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
