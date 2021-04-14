#' @name plot_hist
#'
#' @title This function produces plots for all containing tags
#'
#' @param SangeR object from allign function.
#' @param POI list of points of interesst
#'
#' @return Objekt with plots for the given tags
#'
#' @export
#'


#global Variables
globalVariables(c("position","Values","Samples"))

plot_hist <- function(SangeR, POI = NULL){

  #chromatogramm

  #get values for chromatogram
  #order of channels
  channel_order <- SangeR$abif@data$FWO_.1

  #build data frome for histogramm
  #pool channels in data frame
  data <- data.frame(1:length(SangeR$abif@data$DATA.9), SangeR$abif@data$DATA.9, SangeR$abif@data$DATA.10, SangeR$abif@data$DATA.11, SangeR$abif@data$DATA.12)
  colnames(data) <- c("position",stringr::str_split(channel_order,"")[[1]])

  #add basecalls
  basecalls <- data.frame(SangeR$abif@data$PLOC.1, unlist(strsplit(SangeR$abif@data$PBAS.1, ""))[1:nchar(SangeR$abif@data$PBAS.1)-1])
  colnames(basecalls) <- c("position", "basecall")
  data[basecalls$position,"basecall"] <- basecalls$basecall

  #points of interest

  if (!is.null(POI)){

    POI <- utils::read.table(POI)

    for(point in unlist(POI)){
      chr <- as.numeric(strsplit(point,":")[[1]][1])
      pos <- as.numeric(strsplit(point,":")[[1]][2])
      end <- ifelse(SangeR$ref_pos$strand==-1,SangeR$ref_pos$end_position+SangeR$upstream,SangeR$ref_pos$end_position)
      start <- ifelse(SangeR$ref_pos$strand==1,SangeR$ref_pos$start_position-SangeR$upstream,SangeR$ref_pos$start_position)
        if(chr == SangeR$ref_pos$chromosome_name && pos < end && pos > SangeR $ref_pos$start_position && !any(as.numeric(substr(point,nchar(point)-2,nchar(point))) == substr(SangeR$tags,2,4))){
          tag <- paste0(substr(point,nchar(point)-2,nchar(point)),"wt")
          SangeR$tags <- c(SangeR$tags, tag)
          SangeR$mutations_ref <- c(SangeR$mutations_ref, end - pos)
        }
    }
  }

  #region of interest

  if(length(SangeR$mutations_ref) != 0L){

    #create counter and create vat for used file paths

    PNG_list <- list()
    cnt <- 1

    for(mut in SangeR$mutations_ref){

      #mutation position in abifile
      pos <- mut - SangeR$abi_align@subject@range@start + 1

      #check if position is on the borders of the sequenz
      if(pos > 6 && pos < (length(basecalls$position)-5)){

        #Pick region of interest
        ROI <- data[(SangeR$abif@data$PLOC.1[(pos-5)]):SangeR$abif@data$PLOC.1[(pos+5)],]
        ROI_melt <- reshape2::melt(ROI, id.vars=c("position", "basecall"),  variable.name = "Samples", value.name="Values", na.rm = FALSE)

        #print plot
        plot <- ggplot2::ggplot(ROI_melt, ggplot2::aes(position, Values)) +
                  ggplot2::scale_x_continuous(breaks = data$position[!is.na(data$basecall)], labels =  data$basecall[!is.na(data$basecall)]) +
                  ggplot2::geom_line(ggplot2::aes(color=as.factor(Samples))) +
                  ggplot2::scale_color_discrete(name = "nucleotide") +
                  ggplot2::labs(title = paste0(SangeR$genename, " ", SangeR$tags[cnt])) +
                  ggplot2::theme_bw() +
                  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,))

        #write plot to file
        #grDevices::png(file = paste0("Chromatogramm_", SangeR$Bnummer, "_", SangeR$genename,"_",SangeR$tags[cnt],".png"), width = 1200)
        #print(plot)
        #dev.off()
        PNG_list[[cnt]] <- plot

      } else{
        print("no mutations")
      }
      cnt <- cnt+1
    }
  }  else{

    print("no mutations")

  }
  return(PNG_list)
}
