#' @name run_routine
#'
#' @title This function produces two plots when started in nextflow
#'
#' @param Files
#' @param POI
#'
#' @return two pngs and a csv file
#'
#' @export
#'

library(sangeR)
library(gridExtra)

run_routine <- function(Files, POI){

SangeR <- plot_hist(allign(get_ref(read.ab1(filename = Files))),POI = POI)

if(length(SangeR$PNG_list) > 0){
  ggplot2::ggsave(plot = do.call("grid.arrange", c(SangeR$PNG_list, ncol=1)), paste0(SangeR$Bnummer,".png"), width = 350, height = 300, units='mm')
}
write.table(file = paste0(SangeR$Bnummer,SangeR$genename,".csv"), x = "no mutation found",row.names = FALSE, col.names = FALSE, quote = FALSE)
}

args = commandArgs(trailingOnly=TRUE)

file = args[1]

if (length(args) == 1){
run_routine(Files = args[1])
} else {
  run_routine(Files = args[1], POI = args[2])
}
