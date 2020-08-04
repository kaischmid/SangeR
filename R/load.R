#' read.ab1 module
#'
#' This function reads in a given abi1-file
#'
#' @param filename location of a ab1 file.
#' @param
#'
#' @return Objekt with all information taken from the given file.
#'
#' @export



read.ab1 <- function(filename,cutoff = 0.05, min_seq_len = 20, offset = 33){

  #read ab1file with sangerseqR
  file <- sangerseqR::read.abif(filename = filename)

  #extract B-number and genename
  ids <- file@data$SMPL.1

  #extract B-nummer

  Bnummer <- strsplit(ids, "_")[[1]][1]

  #extrat genename

  genename <- strsplit(ids, "_")[[1]][2]
  genename <- gsub('[[:lower:]]', '', genename)

  #tranfer to fasta

  file_name <- paste0(Bnummer, "_", genename, ".fastq")
  invisible(CrispRVariants::abifToFastq(seqname = file_name, fname = file_name, outfname = file_name, cutoff = cutoff, min_seq_len = min_seq_len, offset = offset))

  #load fastq

  fastq <- Biostrings::readDNAStringSet(filename, format = "fastq")

  #return
  return(fastq)

}
