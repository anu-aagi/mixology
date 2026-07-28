#' Homogeneous variances and zero correlation
homo <- function(sigma = 1){
  function(N){
    Matrix::Diagonal(N, sigma^2)
  }
}

#' Heterogeneous variances and zero correlation
hetero <- function(alpha = 2, beta = 1, sigma = NULL){
  function(N){
    if(is.null(sigma)){
      sigma <- 1/rgamma(N, alpha, beta)
    }
    if(!is.vector(sigma) || !is.numeric(sigma) || length(sigma) != N){
      stop("Invalid 'sigma' specification: Not a numeric vector of length N = ", N)
    }
    Matrix::Diagonal(x = sigma)
  }
}

#' Compound symmetry
cs <- function(sigma = 1, r = 0.5){
  function(N){
    G <- Matrix::Diagonal(N) + r
    Matrix::diag(G) <- Matrix::diag(G) - r
    G * sigma^2
  }
}

#' Collect variance functions
#' @export
vcov <- list(
  homo = homo,
  hetero = hetero,
  cs = cs
)
