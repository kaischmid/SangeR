#if statments for heterocygote decision
#' @name heterozygote
#'
#' @title This function checks for heterozygote Nucleotide and gives back the mutation
#'
#' @param hetero heterozygote base R, Y, S, W, K, or M
#' @param nucleotide Nucleotide from reference
#'
#' @return Nucleotide for mutation
#'
#' @export
#'

heterozygote <- function(hetero, nucleotide){
  if (toupper(hetero) == "R"){
    if (toupper(nucleotide) == "A") {
      return("G")
    }
    else{return("A")}
  }
  if (toupper(hetero) == "Y"){
    if (toupper(nucleotide) == "T") {
      return("C")
    }
    else{return("T")}
  }
  if (toupper(hetero) == "S"){
    if (toupper(nucleotide) == "G") {
      return("C")
    }
    else{return("G")}
  }
  if (toupper(hetero) == "W"){
    if (toupper(nucleotide) == "A") {
      return("T")
    }
    else{return("A")}
  }
  if (toupper(hetero) == "K"){
    if (toupper(nucleotide) == "T") {
      return("G")
    }
    else{return("T")}
  }
  if (toupper(hetero) == "M"){
    if (toupper(nucleotide) == "A") {
      return("C")
    }
    else{return("A")}
  }
}
