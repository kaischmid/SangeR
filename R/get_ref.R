#' get_ref module
#'
#' This function loads the reference for a given SangeR object
#'
#' @param SangeR object from read.ab1() function.
#' @param upstream region before gene to be loaded into ref_seq
#' @param host Host to connect to. Defaults to grch37.ensembl.or
#' @param dataset Dataset you want to use. To see the different datasets available within a biomaRt you can e.g. do: mart = useMart('ensembl'), followed by listDatasets(mart). Default: hsapiens_gene_ensembl
#' @param biomart BioMart database name you want to connect to. Possible database names can be retrieved with the function listMarts. Default: ensembl
#'
#' @return Objekt with all information of the reference.
#'
#' @export


#global Variables
globalVariables(c("exit"))

get_ref <- function(SangeR, upstream = 500, host = "grch37.ensembl.org", dataset = "hsapiens_gene_ensembl", biomart = "ensembl"){

  #define mart
  mart = biomaRt::useMart(biomart = biomart, dataset = dataset, host = host)

  #write upstream to SangeR object
  SamgeR$upstream <- upstream

  #get reference sequence
  SangeR$ref_seq <- ref_sequence <- biomaRt::getSequence(id = SangeR$genename, mart = mart, type = "hgnc_symbol", seqType = "gene_exon_intron", upstream = upstream)

  #ref aminoacid sequence
  temp_aminoacid <- biomaRt::getSequence(id = SangeR$genename, mart = mart, type = "hgnc_symbol", seqType = c("peptide","refseq_peptide",'end_position','start_position'))
  temp_aminoacid$refseq_peptide[temp_aminoacid$refseq_peptide==""] <- NA
  SangeR$ref_amino <- ref_aminoacid <- temp_aminoacid[!is.na(temp_aminoacid$refseq_peptide),]

  #reference position
  SangeR$ref_pos <- ref_position <- biomaRt::getBM(c('hgnc_symbol','chromosome_name','start_position','end_position','strand'),"hgnc_symbol",SangeR$genename, mart)

  #pepitde informations
  pep_info <- biomaRt::getBM(values = strsplit(SangeR$ref_amino$refseq_peptide,split = ';')[[1]][1], "refseq_peptide", attributes = c('chromosome_name','exon_chrom_start','exon_chrom_end'), mart = mart)
  SangeR$pep_info <- pep_info[order(pep_info$exon_chrom_start),]
  SangeR$pep_info$length <- (SangeR$pep_info$exon_chrom_end - SangeR$pep_info$exon_chrom_start - 3)

  #controll if gene could be found

  if(length(ref_sequence$gene_exon_intron) == 0) {
    print("Gene could not been found, please verify that the gene name is correct")
    exit()
  }
  SangeR$upstream <- upstream

  #return
  return(SangeR)

}
