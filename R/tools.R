#' Calculate effective sample size
#' @export
ess <- function(obj, m = NULL, x = NULL, S = NULL, tol = 1e-10) {
  # Check inputs
  if (!inherits(obj, "cocktail")) {
    stop("`obj` must inherit from class 'cocktail'.")
  }

  if (is.null(obj$data$y)) {
    stop("No trait was detected. Run `add_trait` first.")
  }

  N <- obj$data$N
  M <- obj$data$M
  V <- obj$covariance$V

  # Check subsample
  if (!is.null(S)) {
    if (!is.vector(S) || !is.numeric(S)) {
      stop("Invalid S`: Must be a numeric vector of sample indices.")
    }
    S <- unique(as.integer(S))
    if (any(S < 1 | S > N)) {
      stop("Invalid `S`: All indices must be in 1:N (N = ", N, ").")
    }
    N <- length(S)
    V <- V[S, S]
  } else {
    S <- seq_len(N)
  }
  X <- obj$data$X[S, ]

  # Inverse of VCOV matrix (for correlated noises)
  Vinv <- solve(V)

  # Inverse of diagonal VCOV matrix (for uncorrelated noises)
  Dinv <- Matrix::Diagonal(x = 1 / diag(V))

  # Centering and wittering operator: C(A) = A - (A1)(A1)^T / (1^T A 1)
  C <- function(A, N) {
    ones <- rep(1, N)
    A1 <- A %*% ones
    A - (A1 %*% t(A1)) / drop(t(ones) %*% A1)
  }
  XTCVinvX <- t(X) %*% C(Vinv, N) %*% X
  XTCDinvX <- t(X) %*% C(Dinv, N) %*% X

  # If no x or m is supplied, return XTCVinvX and XTCDinvX
  if (is.null(m) & is.null(x)) {
    stop("Neither `m` nor `x` was be provided,")
  }

  # Calculate per-marker ESS / ESS of user defined contrast (x)
  if (!is.null(m)) {
    if (!is.vector(m) || !is.numeric(m)) {
      stop("Invalid m`: Must be a numeric vector of marker indices.")
    }
    m <- unique(as.integer(m))
    if (any(m < 1 | m > M)) {
      stop("Invalid `m`: All indices must be in 1:M (M = ", M, ").")
    }

    gamma <- diag(XTCVinvX) / diag(XTCDinvX)
    gamma <- gamma[m]
    names(gamma) <- obj$data$markers[m]
  } else if (!is.null(x)) {
    if (!is.numeric(x) && !inherits(x, "Matrix")) {
      stop("Invalid `x`: Must be numeric.")
    }

    if (is.vector(x)) {
      if (length(x) != M) {
        stop("Invalid `x`: Must have length M = ", M, ".")
      }
    } else if (inherits(x, "Matrix") || inherits(x, "matrix")) {
      if (nrow(x) != M) {
        stop("Invalid `x`: Must have row dimension equals to M = ", M, ".")
      }
    } else {
      stop("Invalid `x`: Must be a numeric vector or matrix.")
    }
    x <- Matrix::Matrix(x, nrow = M)

    gamma <- apply(x, 2, function(xj) {
      drop(t(x) %*% XTCVinvX %*% x) /
        drop(t(x) %*% XTCDinvX %*% x)
    })
  }

  gamma
}

#' Calculate or extract the covariance matrix
#' @export
covariance <- function(obj, which = "V", calculate = FALSE) {
  # Check inputs
  if (!inherits(obj, "cocktail")) {
    stop("`obj` must inherit from class 'cocktail'.")
  }

  if (is.null(obj$data$y)) {
    stop("No trait was detected. Run `add_trait` first.")
  }

  N <- obj$data$N
  M <- obj$data$M
  random_effects <- sapply(obj$covariance$G, is.null) %>%
    isFALSE() %>%
    names(.)[.]

  if (which == "V") {
    if (calculate) {
      V <- obj$covariance$R
      if (length(random_effects) > 0) {
        for (g in random_effects) {
          Z <- obj$design$random[[g]]
          G <- obj$covariance$G[[g]]
          V <- V + Z %*% G %*% t(Z)
        }
      }
    } else {
      V <- obj$covariance$V
    }
    return(V)
  }

  if(which == "R") {
    return(objcovariance$R)
  }

  if (which %in% random_effects) {
    return(objcovariance$G[[which]])
  } else {
    stop("Invalid `which` specification: Unknown variance component.")
  }

}

#' Extract genotype data
#' @export
genotype <- function(obj, n = NULL, m = NULL) {
  # Check inputs
  if (!inherits(obj, "cocktail")) {
    stop("`obj` must inherit from class 'cocktail'.")
  }

  N <- obj$data$N
  M <- obj$data$M

  if (!is.null(m)) {
    if (!is.vector(m) || !is.numeric(m)) {
      stop("Invalid m`: Must be a numeric vector of marker indices.")
    }
    m <- unique(as.integer(m))
    if (any(m < 1 | m > M)) {
      stop("Invalid `m`: All indices must be in 1:M (M = ", M, ").")
    }
  } else {
    m <- seq_len(M)
  }

  if (!is.null(n)) {
    if (!is.vector(n) || !is.numeric(n)) {
      stop("Invalid n`: Must be a numeric vector of marker indices.")
    }
    n <- unique(as.integer(n))
    if (any(n < 1 | n > M)) {
      stop("Invalid `n`: All indices must be in 1:N (N = ", N, ").")
    }
  } else {
    n <- seq_len(N)
  }


  obj$data$X[n,M]
}

#' Extract metadata
#' @export
metadata <- function(obj, which = NULL) {
  # Check inputs
  if (!inherits(obj, "cocktail")) {
    stop("`obj` must inherit from class 'cocktail'.")
  }

  metas <- names(obj$data$meta)

  if (is.null(which)) which <- metas

  obj$data$meta[which]
}

