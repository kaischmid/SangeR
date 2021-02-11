#' @name translate
#'
#' @title This function takes the Aa exchange and reference to give back the Aa mutation
#'
#' @param Aa AminoAcid from the reference
#' @param exchange result from heterozygote from the mutation position
#' @param reference the reference with the two position before and after the mutation
#'
#' @return Objekt with plots for the given tags
#'
#' @export
#'

translate <- function(Aa, exchange, reference){

  if (seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 0, sens = "F") == Aa){
    substr(reference, 3, 3) <- exchange
    return(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 0, sens = "F"))
  }
  if(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 1, sens = "F") == Aa){
    substr(reference, 3, 3) <- exchange
    return(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 1, sens = "F"))
  }
  if(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 2, sens = "F") == Aa){
    substr(reference, 3, 3) <- exchange
    return(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 2, sens = "F"))
  }
  if(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 0, "R") == Aa){
    substr(reference, 3, 3) <- exchange
    return(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 0, sens = "R"))
  }
  if(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 1, sens = "R") == Aa){
    substr(reference, 3, 3) <- exchange
    return(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 1, sens = "R"))
  }
  if(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 2, sens = "R") == Aa){
    substr(reference, 3, 3) <- exchange
    return(seqinr::translate(as.vector(tolower(strsplit(reference,"")[[1]])), frame = 2, sens = "R"))
  }
}
