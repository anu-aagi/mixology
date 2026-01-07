#' No transform
ident <- function(){
  function(obj){
    obj$data$X
  }
}

#' VanRaden transform
vanraden <- function(){
  function(obj){
    W <- obj$data$X
    N <- obj$data$N
    pm <- colSums(W) / (2 * N) # allele freq per marker (diploid X)
    pm <- pmin(pmax(pm, 1e-6), 1 - 1e-6) # guard against 0 or 1

    W <- sweep(W, 2, 2 * pm, "-")
    W <- sweep(W, 2, sqrt(2 * pm * (1 - pm)), "/")
    W
  }
}

#' Collect transform functions
#' @export
feature_map <- list(
  ident = ident,
  vanraden = vanraden
)
