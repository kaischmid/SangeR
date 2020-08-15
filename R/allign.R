#' allign module
#'
#' This function loads the reference for a given SangeR object
#'
#' @param SangeR object from get_ref() function.
#' @param upstream region before gene to be loaded into ref_seq
#' @param host Host to connect to. Defaults to grch37.ensembl.or
#' @param dataset Dataset you want to use. To see the different datasets available within a biomaRt you can e.g. do: mart = useMart('ensembl'), followed by listDatasets(mart). Default: hsapiens_gene_ensembl
#' @param biomart BioMart database name you want to connect to. Possible database names can be retrieved with the function listMarts. Default: ensembl
#'
#' @return Objekt with all information about the alignment.
#'
#' @export

allign <- function(SangeR){

  #load sequence in variable
  sequence <- toString(SangeR$fastq[1])

  #allignments

  mart_align_forward <- Biostrings::pairwiseAlignment(subject = Biostrings::DNAString(SangeR$ref_seq$gene_exon_intron), pattern = Biostrings::DNAString(sequence), type = "global-local")
  abi_align_forward <- Biostrings::pairwiseAlignment(subject = Biostrings::DNAString(SangeR$ref_seq$gene_exon_intron), pattern = Biostrings::DNAString(SangeR$abif@data$PBAS.1), type = "global-local")
  mart_align_reverse <- Biostrings::pairwiseAlignment(subject = Biostrings::reverseComplement(Biostrings::DNAString(SangeR$ref_seq$gene_exon_intron)), pattern = Biostrings::DNAString(sequence), type = "global-local")
  abi_align_reverse <- Biostrings::pairwiseAlignment(subject = Biostrings::reverseComplement(Biostrings::DNAString(SangeR$ref_seq$gene_exon_intron)), pattern = SangeR$abif@data$PBAS.1, type = "global-local")

  #choose if sequence is on forward or reverse strand

  if(mart_align_forward@score>=mart_align_reverse@score){

    abi_align <- abi_align_forward
    mart_align <- mart_align_forward
    SangeR$strand <- "forward"

  } else {

    abi_align <- abi_align_reverse
    mart_align <- mart_align_reverse
    SangeR$strand <- "reverse"

  }

  #Mutation

  #write tags for mismatches

  mutations <- c()
  align <- c()

  #use only common shared mismatches between abi file and the fasta base call

  mutations <- mart_align@subject@mismatch[[1]][mart_align@subject@mismatch[[1]] %in% abi_align@subject@mismatch[[1]]]
  align <- mart_align

  #combine tag elements for each string

  if(length(align@subject@mismatch[[1]]) != 0){

    #write the tags

    cnt <- 1
    tags <- c()

    for(mut in mutations[[1]]){

      #find positions

      if(SangeR$ref_pos$strand == -1) {

        #mutation position for - strand
        mutpos <- paste0("chr", SangeR$ref_pos$chromosome_name,":", SangeR$ref_pos$end_position - mut + SangeR$upstream + 1)

      } else {

        #mutation position for + strand
        mutpos <- paste0("chr",SangeR$ref_pos$chromosome_name,":", SangeR$ref_pos$start_position + mut + SangeR$upstream + 1)
      }
      chr_pos <- strsplit(mutpos,":")[[1]][2]

      #check if mutation is in CDS
      if(any((SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos))) {

        #write tag for protein region
        if(SangeR$ref_pos$strand == -1){

          #for reverse strand

          pos <- (((as.numeric(Sanger$pep_info$exon_chrom_end[(SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos)])-as.numeric(chr_pos)-3) + sum(SangeR$pep_info$length[(which(((SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos)))+1):length(SangeR$pep_info$length)]))/3)

          Aa <- stringr::str_sub(ref_aminoacid$peptide, align@subject@mismatch[[1]][cnt])
          Aa_mut <- stringr::str_sub(sequence, align@pattern@mismatch[[1]][cnt], align@pattern@mismatch[[1]][cnt])
          tags <- c(tags, paste0(Aa, Aa_mut))

        }else{

          #for forward strand

          Aa <- stringr::str_sub(ref_aminoacid$peptide, align@subject@mismatch[[1]][cnt])
          Aa_mut <- stringr::str_sub(sequence, align@pattern@mismatch[[1]][cnt], align@pattern@mismatch[[1]][cnt])
          tags <- c(tags, paste0(Aa, Aa_mut))
        }

      } else {

        #write tag for off-region
        base <- stringr::str_sub(SangeR$ref_seq$gene_exon_intron, mut, mut)
        mut <- stringr::str_sub(sequence, align@pattern@mismatch[[1]][cnt], align@pattern@mismatch[[1]][cnt])
        tags <- c(tags, paste0(base,stringr::str_sub(mutpos, -3, -1),mut))

      }
      #for counter
      cnt <- cnt+1
    }
    SangeR$tags <- tags
  }

  return(SangeR)

}
