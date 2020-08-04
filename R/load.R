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



read.ab1 <- function(filename){

  #read ab1file with sangerseqR
  file <- sangerseqR::read.abif(filename = abi_file)

  #extract B-number and genename
  ids <- file@data$SMPL.1

  #extract B-nummer

  Bnummer <- strsplit(ids, "_")[[1]][1]

  #extrat genename

  genename <- strsplit(ids, "_")[[1]][2]
  genename <- gsub('[[:lower:]]', '', genename)

  #tranfer to fasta

  filename <- paste0(Bnummer, "_", genename, ".fastq")
  invisible(CrispRVariants::abifToFastq(seqname = filename, fname = abi_file, outfname = filename, cutoff = cutoff, min_seq_len = min_seq_len, offset = offset))

  #return
  return(objects)

}
