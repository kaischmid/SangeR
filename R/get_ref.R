#' @name get_ref
#'
#' @title This function loads the reference for a given SangeR object
#'
#' @param SangeR object from read.ab1() function.
#' @param upstream region before gene to be loaded into ref_seq
#' @param host Host to connect to. Defaults to grch37.ensembl.or
#' @param dataset Dataset you want to use. To see the different datasets available within a biomaRt you can e.g. do: mart = useMart("ensembl"), followed by listDatasets(mart). Default: hsapiens_gene_ensembl
#' @param biomart BioMart database name you want to connect to. Possible database names can be retrieved with the function listMarts. Default: ensembl
#'
#' @return Objekt with all information of the reference.
#'
#' @export

get_ref <- function(SangeR, upstream = 500, host = "grch37.ensembl.org", dataset = "hsapiens_gene_ensembl", biomart = "ensembl"){

  #define mart
  mart = biomaRt::useMart(biomart = biomart, dataset = dataset, host = host)

  #write upstream to SangeR object
  SangeR$upstream <- upstream

  #get reference sequence
  try(SangeR$ref_seq <- biomaRt::getSequence(id = SangeR$genename, mart = mart, type = "hgnc_symbol", seqType = "gene_exon_intron", upstream = upstream), silent = TRUE)

  #control if gene could be found

  if(length(SangeR$ref_seq$gene_exon_intron) == 0){
    cnt <- 1
    SangeR$Bnummer <- SangeR$Bnummer
    while(length(SangeR$ref_seq$gene_exon_intron) == 0 & cnt < 5) {
      genename <- strsplit(SangeR$ids, SangeR$delimiter)[[1]][cnt]
      genename <- gsub('[[:lower:]]', '', genename)
      genename <- gsub('-','',genename)
      SangeR$genename <- genename
      try(SangeR$ref_seq <- biomaRt::getSequence(id = SangeR$genename, mart = mart, type = "hgnc_symbol", seqType = "gene_exon_intron", upstream = upstream), silent = TRUE)
      cnt <- cnt +1
    }
  }

  if(length(SangeR$ref_seq$gene_exon_intron) == 0) {
    stop(paste0("Gene could not been found, please verify that the gene name is correct\nYou entered the genename: ",SangeR$genename))
  }

  #ref aminoacid sequence
  temp_aminoacid <- biomaRt::getSequence(id = SangeR$genename, mart = mart, type = "hgnc_symbol", seqType = c("peptide","refseq_peptide","end_position","start_position"))
  temp_aminoacid$refseq_peptide[temp_aminoacid$refseq_peptide==""] <- NA
  SangeR$ref_amino <- temp_aminoacid[!is.na(temp_aminoacid$refseq_peptide),]

  #reference position
  SangeR$ref_pos <- biomaRt::getBM(c("hgnc_symbol","chromosome_name","start_position","end_position","strand"), "hgnc_symbol", SangeR$genename, mart)

  #pepitde informations

  pep_info <- biomaRt::getBM(values = strsplit(SangeR$ref_amino$refseq_peptide,split = ";")[[1]][1], "refseq_peptide", attributes = c("chromosome_name","cdna_coding_start","cdna_coding_end","exon_chrom_start","exon_chrom_end"), mart = mart)
  SangeR$pep_info <- pep_info[order(pep_info$cdna_coding_start),]
  SangeR$pep_info$length <- (SangeR$pep_info$cdna_coding_end - SangeR$pep_info$cdna_coding_start) + 1


  #return
  return(SangeR)

}
