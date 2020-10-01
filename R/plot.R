#' plot_hist module
#'
#' This function produces plots for all containing tags
#'
#' @param SangeR object from allign function.
#'
#' @return Objekt with plots for the given tags
#'
#' @export
#'

#global Variables
globalVariables(c("mutations","abi_align","mart_align","position","Values","Samples","genename","tags"))

plot_hist <- function(SangeR){
#chromatogramm

  #get values for chromatogram
  #order of channels
  channel_order <- SangeR$abi@data$FWO_.1

  #build data frome for histogramm
  #pool channels in data frame
  data <- data.frame(1:length(SangeR$abi@data$DATA.9), SangeR$abi@data$DATA.9, SangeR$abi@data$DATA.10, SangeR$abi@data$DATA.11, SangeR$abi@data$DATA.12)
  colnames(data) <- c("position",stringr::str_split(channel_order,"")[[1]])

  #add basecalls
  basecalls <- data.frame(SangeR$abi@data$PLOC.1, unlist(strsplit(SangeR$abi@data$PBAS.1, ""))[1:nchar(SangeR$abi@data$PBAS.1)-1])
  colnames(basecalls) <- c("position", "basecall")
  data[basecalls$position,"basecall"] <- basecalls$basecall


  #region of interest

  if(is.integer(mutations[[1]]) && length(mutations[[1]]) != 0L){

    #create counter and create vat for used file paths

    PNG_list <- c()
    cnt <- 1

    for(mut in mutations[[1]]){

      #mutation position in abifile
      pos <- as.numeric(abi_align@pattern@mismatch[which(mart_align@subject@mismatch %in% abi_align@subject@mismatch)]) #error in pick

      #check if position is on the borders of the sequenz
      if(pos > 6 && pos < (length(basecalls$position)-5)){

        #Pick region of interest
        ROI <- data[(SangeR$abi@data$PLOC.1[(pos-5)]):SangeR$abi@data$PLOC.1[(pos+5)],]
        ROI_melt <- reshape2::melt(ROI, id.vars=c("position", "basecall"),  variable.name = "Samples", value.name="Values", na.rm = FALSE)

        #print plot
        plot <- ggplot2::ggplot(ROI_melt, ggplot2::aes(position, Values)) +
                  ggplot2::scale_x_continuous(breaks = data$position[!is.na(data$basecall)], labels =  data$basecall[!is.na(data$basecall)]) +
                  ggplot2::geom_line(ggplot2::aes(color=as.factor(Samples))) +
                  ggplot2::scale_color_discrete(name = "nucleotide") +
                  ggplot2::labs(title = paste0(genename, " ", tags)) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,))

        #write plot to file
        grDevices::png(file = paste0("Chromatogramm_", SangeR$Bnummer, "_", SangeR$genename,".png"), width = 1200)
        print(plot)
        dev.off()
        PNG_list <- c(PNG_list, paste0(getwd(),"/Chromatogramm_", SangeR$Bnummer, "_", SangeR$genename,".png"))
      } else{
        print("no mutations")
      }
      cnt <- cnt+1
    }
  }  else{

    print("no mutations")

  }
}
