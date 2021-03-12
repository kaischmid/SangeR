
library(biomaRt)
library(Homo.sapiens)
library(GenomicRanges)
host = "grch37.ensembl.org"
dataset = "hsapiens_gene_ensembl"
biomart = "ensembl"


#read parameters
args = commandArgs(trailingOnly=TRUE)

#read in
bed <- read.table("../../Downloads/S31285117_Covered.bed", sep = "\t",skip = 2, stringsAsFactors=FALSE, quote="")

#biomart

mart = biomaRt::useMart(biomart = biomart, dataset = dataset, host = host)

region <- lapply(1:length(bed$V1), function(i){
  list(substr(bed[i,]$V1,1,4),bed[i,]$V2,bed[i,]$V3)
})

result_hgnc <- biomaRt::getBM(values = unlist(region[1:500]), useCache = FALSE, mart = mart, attributes = c("hgnc_symbol"), filters = "chromosomal_region",verbose = TRUE, bmHeader = TRUE, )

result_hgnc_2 <- biomaRt::getBM(values = paste0(substr(bed[i,]$V1,4,4),":",bed[i,]$V2,",",substr(bed[i,]$V2,4,4),":",bed[i,]$V3), mart = mart, attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"))


#granges attempt

#read in
bed <- read.table("../../Downloads/S31285117_Covered.bed", sep = "\t",skip = 2, stringsAsFactors=FALSE, quote="")

region <- lapply(1:length(bed$V1), function(i){
  list(substr(bed[i,]$V1,1,4),bed[i,]$V2,bed[i,]$V3)
})

geneRanges <- function(db, column="SYMBOL")
{
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), lengths(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

gns <- suppressMessages(geneRanges(Homo.sapiens))

splitColumnByOverlap <-  function(query, subject, column="SYMBOL", ...)
{
  olaps <- findOverlaps(query, subject, ...)

  f1 <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))

  splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}



table <- data.frame(t(as.data.frame(matrix(unlist(region),nrow = 3))))

colnames(table) <- c("chr", "start", "end")

cnv <- makeGRangesFromDataFrame(table)


gene.Symbol <- suppressWarnings(splitColumnByOverlap(gns, cnv, "SYMBOL"))

symbol_temp <- as.vector(gene.Symbol)

symbol <- vapply(symbol_temp, paste, collapse = ", ", character(1L))


table$symbol <- symbol

write.table(table, file = "hgnc.tsv", sep = "\t", row.names = FALSE)
write.table(unlist(result)[1:500], file = "bed.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

