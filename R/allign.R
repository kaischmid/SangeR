#' @name  allign
#'
#' @title  This function loads the reference for a given SangeR object
#'
#' @param SangeR object from get_ref() function.
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

    SangeR$abi_align <- abi_align_forward
    SangeR$mart_align <- mart_align_forward
    SangeR$align_seq <- Biostrings::DNAString(SangeR$ref_seq$gene_exon_intron)
    SangeR$strand <- "forward"

  } else {

    SangeR$abi_align <- abi_align_reverse
    SangeR$mart_align <- mart_align_reverse
    SangeR$align_seq <- Biostrings::reverseComplement(Biostrings::DNAString(SangeR$ref_seq$gene_exon_intron))
    SangeR$strand <- "reverse"

  }

  #Mutation

  #use only common shared mismatches between abi file and the fasta base call

  SangeR$mutations_ref <- SangeR$mart_align@subject@mismatch[[1]][SangeR$mart_align@subject@mismatch[[1]] %in% SangeR$abi_align@subject@mismatch[[1]]]
  if(SangeR$param != "default"){SangeR$mutations_ref <- SangeR$mart_align@subject@mismatch[[1]]}
  SangeR$mutations_abi <- SangeR$mart_align@pattern@mismatch[[1]][SangeR$mart_align@subject@mismatch[[1]] %in% SangeR$abi_align@subject@mismatch[[1]]]
  if(SangeR$param != "default"){SangeR$mutations_abi <- SangeR$mart_align@pattern@mismatch[[1]]}
  SangeR$align <- SangeR$mart_align

  #combine tag elements for each string

  if(length(SangeR$mutations_ref) != 0){

    #write the tags

    cnt <- 1
    tags <- c()

      for(mut in SangeR$mutations_ref[[1]]){

      #find positions

      if(SangeR$strand == "forward") {

        #mutation position for forward strand
        mutpos <- paste0("chr", SangeR$ref_pos$chromosome_name,":", SangeR$ref_pos$end_position - mut + SangeR$upstream + 1)

      } else {

        #mutation position for reward strand
        mutpos <- paste0("chr",SangeR$ref_pos$chromosome_name,":", SangeR$ref_pos$start_position + mut - 1)
      }
      chr_pos <- strsplit(mutpos,":")[[1]][2]

      #check if mutation is in CDS
      if(any((SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos))) {

        #write tag for protein region
        if(SangeR$strand == "forward"){

          #for forward strand

          boolean <- (SangeR$pep_info$exon_chrom_end<chr_pos)
          boolean[which((SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos))] <- FALSE

          AA_pos <- ceiling((sum(SangeR$pep_info$length[boolean]) + (as.numeric(chr_pos) - SangeR$pep_info$exon_chrom_start[(SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos)]))/3)

          exchange <- if(!substr(SangeR$fastq, SangeR$mutations_abi, SangeR$mutations_abi) %in% c("A","C","T","G")){heterozygote(substr(SangeR$fastq, SangeR$mutations_abi, SangeR$mutations_abi), stringr::str_sub(SangeR$align_seq, mut, mut))
            }else {substr(SangeR$fastq, SangeR$mutations_abi, SangeR$mutations_abi)}
          Aa_mut <- translate(Aa, exchange, stringr::str_sub(SangeR$align_seq, mut-2, mut+2))

          tags <- c(tags, paste0(Aa, sprintf("%03d", AA_pos), Aa_mut))

        }else{

          #for reverse strand

          boolean <- (SangeR$pep_info$exon_chrom_end>chr_pos)
          boolean[which((SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos))] <- FALSE

          AA_pos <- ceiling((sum(SangeR$pep_info$length[boolean]) + (SangeR$pep_info$exon_chrom_end[(SangeR$pep_info$exon_chrom_start<chr_pos) == (SangeR$pep_info$exon_chrom_end>chr_pos)]) - as.numeric(chr_pos))/3)

          Aa <- stringr::str_sub(SangeR$ref_amino$peptide, AA_pos, AA_pos)

          exchange <- if(!substr(SangeR$fastq, SangeR$mutations_abi, SangeR$mutations_abi) %in% c("A","C","T","G")){heterozygote(substr(SangeR$fastq, SangeR$mutations_abi, SangeR$mutations_abi), stringr::str_sub(SangeR$align_seq, mut, mut))
            }else {substr(SangeR$fastq, SangeR$mutations_abi, SangeR$mutations_abi)}

          Aa_mut <- translate(Aa, exchange, stringr::str_sub(SangeR$align_seq, mut-2, mut+2))

          tags <- c(tags, paste0(Aa, sprintf("%03d", AA_pos), Aa_mut))
        }

      } else {

        #write tag for off-region
        base <- stringr::str_sub(SangeR$ref_seq$gene_exon_intron, mut, mut)
        mutation <- stringr::str_sub(SangeR$fastq, SangeR$mutations_abi, SangeR$mutations_abi)
        tags <- c(tags, paste0(base,stringr::str_sub(mutpos, -3, -1),mutation))

      }
      #for counter
      cnt <- cnt+1
    }
    SangeR$tags <- tags
  }

  return(SangeR)

}
