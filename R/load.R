#' @name read.ab1
#'
#' @title This function reads in a given abi1-file.
#'
#' @param filename location of a ab1 file.
#' @param delimiter delimiter in genename/ID of ab1 file
#' @param ID_pos position of ID in genename/ID of ab1 file
#' @param genename gene name of analysed gene in ab1 file
#' @param genename_pos position of genename in genename/ID of ab1 file
#' @param cutoff cutoff for Basecall of the fastq. Default: 0.05
#' @param min_seq_len minimum sequence length for Basecall of the fastq. Default: 20
#' @param offset Offset for Basecall of the fastq. Default: 33
#'
#' @return Objekt with all information genename,bnummer,filename,abif,fastq.
#'
#' @export



read.ab1 <- function(filename, delimiter = "_", ID_pos = 1, genename = "", genename_pos = 2,cutoff = 0.05, min_seq_len = 20, offset = 33){

  #read ab1file with sangerseqR
  abif <- sangerseqR::read.abif(filename = filename)

  #extract B-number and genename
  ids <- abif@data$SMPL.1

  #extract B-nummer

  Bnummer <- strsplit(ids, delimiter)[[1]][ID_pos]

  #extract genename

  if(genename == ""){
    genename <- strsplit(ids, delimiter)[[1]][genename_pos]
    genename <- gsub('[[:lower:]]', '', genename)
    genename <- gsub('-','',genename)
  } else {genename <- genename}


  #tranfer to fasta

  file_name <- paste0(Bnummer, "_", genename, ".fastq")
  invisible(CrispRVariants::abifToFastq(seqname = file_name, trim = TRUE, fname = filename, outfname = file_name, cutoff = cutoff, min_seq_len = min_seq_len, offset = offset))

  #load fastq

  fastq <- Biostrings::readDNAStringSet(file_name, format = "fastq")

  #clean up

  file.remove(file_name)

  #param meter of read in

  if(cutoff == 0.05 && min_seq_len == 20 && offset == 33){param <- NULL
}else {param <- c(cutoff, min_seq_len, offset)}

  object <- function(g,b,f,a,fa,p,i,d,o,m,c) {
    value <- list(genename=g, Bnummer=b, filename=f, abif=a, fastq=fa, param=p, ids=i, delimiter=d,offset=offset,min_seq_len=min_seq_len,cutoff=cutoff)
    attr(value, "class") <- "SangeR"
    value
  }


  #return
  return(object(genename,Bnummer,filename,abif,fastq,param,ids,delimiter,offset,min_seq_len,cutoff))

}
