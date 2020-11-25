#' read.ab1 module
#'
#' This function reads in a given abi1-file
#'
#' @param filename location of a ab1 file.
#' @param cutoff cutoff for Basecall of the fastq. Default: 0.05
#' @param min_seq_len minimum sequence length for Basecall of the fastq. Default: 20
#' @param offset Offset for Basecall of the fastq. Default: 33
#'
#' @return Objekt with all information genename,bnummer,filename,abif,fastq.
#'
#' @export



read.ab1 <- function(filename, cutoff = 0.05, min_seq_len = 20, offset = 33){



  #read ab1file with sangerseqR
  abif <- sangerseqR::read.abif(filename = filename)

  #extract B-number and genename
  ids <- abif@data$SMPL.1

  #extract B-nummer

  Bnummer <- strsplit(ids, "_")[[1]][1]

  #extrat genename

  genename <- strsplit(ids, "_")[[1]][1]
  genename <- gsub('[[:lower:]]', '', genename)

  #tranfer to fasta

  file_name <- paste0(Bnummer, "_", genename, ".fastq")
  invisible(CrispRVariants::abifToFastq(seqname = file_name, trim = TRUE, fname = filename, outfname = file_name, cutoff = cutoff, min_seq_len = min_seq_len, offset = offset))

  #load fastq

  fastq <- Biostrings::readDNAStringSet(file_name, format = "fastq")

  object <- function(g,b,f,a,fa) {
    value <- list(genename=g, Bnummer=b, filename=f, abif=a, fastq=fa)
    attr(value, "class") <- "SangeR"
    value
  }


  #return
  return(object(genename,Bnummer,filename,abif,fastq))

}
