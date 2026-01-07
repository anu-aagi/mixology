#' Add a trait
#' @export
add_trait <- function(obj,
                      beta0 = 0,
                      R = 1,
                      scale = FALSE,
                      seed = obj$seed) {
  # Save call
  obj$params$add_pheno <- .capture_all_args(add_trait, drop = c("obj"))

  # Check inputs
  if (!inherits(obj, "cocktail")) {
    stop("`obj` must inherit from class 'cocktail'.")
  }

  if (is.numeric(seed)) set.seed(seed)

  N <- obj$data$N

  # Noise component
  if (is.vector(R)) {
    if (any(length(R) %in% c(1, N)) && is.numeric(R)) {
      eps <- rnorm(N, 0, sqrt(R))
      R <- Matrix::Diagonal(n = N, x = R)
    } else {
      stop("Invalid `R` specification: A vector deteced, but not a numeric of length N = ", N)
    }
  } else {
    m <- .check.vcov.matrix(N, R)
    if (!is.null(m)) {
      stop("Invalid `R` specification: ", m)
    }

    R <- Matrix::Matrix(R)
    L <- t(chol(R))
    u <- rnorm(N, 0, 1)
    eps <- as.numeric(L %*% u)
  }

  trRC <- sum(diag(R) - rowSums(R) / N)
  expected_var <- trRC / (N - 1)

  # Scale to match the desired variance
  if (scale) {
    eps <- eps * sqrt(expected_var) / sd(eps)
  }

  # save results
  if (is.null(obj$data$y)) {
    message("No trait was detected. Creating a new trait.")

    # Create component matrix
    obj$component <- matrix(0,
                            nrow = N,
                            ncol = 2,
                            dimnames = list(
                              obj$data$samples,
                              c("Intercept", "Noise")
                            )
    ) %>% data.frame()

    # Create summary statistics
    obj$summary$expected_var_composition <-
      obj$summary$sample_var_composition <- c(Noise = 0)
    obj$summary$mean <- c(Expected = 0, Sample = 0)
  } else {
    message("Updating the existing trait.")
  }

  # Save simulations
  obj$covariance$R <- R
  obj$component$Intercept <- rep(beta0, N)
  obj$component$Noise <- eps
  obj$data$y <- rowSums(obj$component)

  # Calculate summary stats
  obj$summary$mean["Expected"] <- beta0
  obj$summary$mean["Sample"] <- mean(obj$data$y)
  obj$summary$expected_var_composition["Noise"] <- expected_var
  obj$summary$sample_var_composition["Noise"] <- var(eps)
  obj$summary$expected_total_var <- sum(obj$summary$expected_var_composition)
  obj$summary$sample_total_var <- var(obj$data$y)

  # Calculate covariance v
  obj$covariance$v <- covariance(obj, calculate = TRUE)

  obj
}

#' Add additive effect
#' @export
add_additive <- function(obj,
                         EQ = 10,
                         r = 1,
                         b = NULL,
                         scale = FALSE,
                         seed = NULL) {
  # Save call
  obj$params$add_additive <- .capture_all_args(add_additive, drop = c("obj"))

  # Check inputs
  if (!inherits(obj, "cocktail")) {
    stop("`obj` must inherit from class 'cocktail'.")
  }

  if (is.null(obj$data$y)) {
    stop("No trait was detected. Run `add_trait` first.")
  }
  if (!is.null(obj$component$Additive)) {
    message("Updating the existing additive effect.")
  }

  if (is.numeric(seed)) set.seed(seed)


  # Genotype data
  X <- obj$data$X
  N <- obj$data$N
  M <- obj$data$M

  # Resolve Laplace rate parameter b
  trXtCX <- sum(X^2) - sum(t(X) * colMeans(X))

  if (!is.null(b)) {
    message("`b` is provided, suppress `r`.")
    expected_var <- 2 * b^2 * EQ * trXtCX / M / (N - 1)
  } else {
    if (!is.numeric(r) || length(r) != 1 || r <= 0) {
      stop("Invalid `r` specification: Must be a single numeric.")
    }
    var_noise <- obj$summary$expected_var_composition["Noise"]
    expected_var <- r * var_noise
    b <- sqrt(expected_var * M * (N - 1) / (2 * EQ * trXtCX))
  }

  # Simulate additive effect
  p <- min(1, EQ / M)
  Q <- rbinom(1, M, p) # causal markers (≤ G)

  # Laplace sampler
  rlaplace <- function(n, b) {
    sample(c(-1, 1), n, replace = TRUE) * rexp(n, rate = 1 / b)
  }

  beta <- numeric(M) # initialize effects
  names(beta) <- obj$data$markers # add names
  S_Q <- sample(1:M, Q) # causal marker indices
  beta[S_Q] <- rlaplace(Q, b) # assign effects (consider put distance restriction
  # on causal markers? i.e., they can't be next to each other)

  # Additive component
  additive <- as.numeric(X %*% beta)

  # Scale to match the desired variance
  if (scale) {
    scale_additive <- sqrt(expected_var) / sd(additive)
    additive <- additive * scale_additive
    # Scale effect vectors
    beta <- beta * scale_additive
  }

  # Save simulations
  obj$data$added["Additive"] <- "fixed"
  obj$design$fixed$Additive <- X
  obj$effects$fixed$Additive <- beta
  obj$component$Additive <- additive
  obj$data$y <- rowSums(obj$component)

  # Calculate summary stats
  obj$summary$mean["Sample"] <- mean(obj$data$y)
  obj$summary$expected_var_composition["Additive"] <- expected_var
  obj$summary$sample_var_composition["Additive"] <- var(additive)
  obj$summary$expected_total_var <- sum(obj$summary$expected_var_composition)
  obj$summary$sample_total_var <- var(obj$data$y)

  obj
}

#' Add polygenic effect
#' @export
add_polygenic <- function(obj,
                          W = feature_map$vanraden(),
                          G_Polygenic = 1,
                          r = NULL,
                          scale = FALSE,
                          seed = NULL) {
  # Save call
  obj$params$add_polygenic <- .capture_all_args(add_polygenic, drop = c("obj"))

  # Check inputs
  if (!inherits(obj, "cocktail")) {
    stop("`obj` must inherit from class 'cocktail'.")
  }

  if (is.null(obj$data$y)) {
    stop("No trait was detected. Run `add_trait` first.")
  }
  if (!is.null(obj$effects$Polygenic)) {
    message("Updating the existing polygenic effect.")
  }

  if (is.numeric(seed)) set.seed(seed)

  N <- obj$data$N

  # Standardized genotype matrix
  if (is.function(W)) {
    W <- W(obj)
  }
  if (!(inherits(W, "Matrix") || inherits(W, "matrix"))) {
    stop("Invalid 'W' specification: Not a matrix or a Matrix object.")
  }
  if (nrow(W) != N) {
    stop("Invalid 'W' specification: Incorrect dimension.")
  }

  rownames(W) <- obj$data$samples
  M <- ncol(W)

  # Variance component
  if (is.vector(G_Polygenic)) {
    if (any(length(G_Polygenic) %in% c(1, M)) && is.numeric(G_Polygenic)) {
      G_Polygenic <- Matrix::Diagonal(n = M, x = G_Polygenic)
    } else {
      stop("Invalid `G_Polygenic` specification: A vector deteced, but not a numeric of length N = ", N)
    }
  } else {
    m <- .check.vcov.matrix(M, G_Polygenic)
    if (!is.null(m)) {
      stop("Invalid `G_Polygenic` specification: ", m)
    }
    G_Polygenic <- Matrix::Matrix(G_Polygenic)
  }

  WG <- W %*% G_Polygenic
  trWVWTC <- sum(rowSums(WG * W) - WG %*% colSums(W) / N)

  # Calculate expected variance
  if (!is.null(r)) {
    if (!is.numeric(r) || length(r) != 1 || r <= 0) {
      stop("Invalid `r` specification: Must be a single numeric.")
    }
    message("`r` is provided, rescale the variance component.")
    var_noise <- obj$summary$expected_var_composition["Noise"]
    expected_var <- r * var_noise
    G_Polygenic <- expected_var * G_Polygenic * (N - 1) / trWVWTC
  } else {
    expected_var <- trWVWTC / (N - 1)
  }

  # Simulate polygenic effect
  L <- t(chol(G_Polygenic))
  u <- rnorm(M, 0, 1)
  u <- as.numeric(L %*% u)
  names(u) <- colnames(W) # Add name

  # Polygenic component
  polygenic <- as.numeric(W %*% u)

  # Scale to match the desired variance
  if (scale) {
    scale_polygenic <- sqrt(expected_var) / sd(polygenic)
    polygenic <- polygenic * scale_polygenic
    # Scale effect vectors
    u <- u * scale_polygenic
  }

  # Save simulations
  obj$data$added["Polygenic"] <- "random"
  obj$effects$random$Polygenic <- u
  obj$component$Polygenic <- polygenic
  obj$data$W <- W
  obj$design$random$Polygenic <- W
  obj$covariance$G$Polygenic <- G_Polygenic
  obj$data$y <- rowSums(obj$component)

  # Calculate summary stats
  obj$summary$mean["Sample"] <- mean(obj$data$y)
  obj$summary$expected_var_composition["Polygenic"] <- expected_var
  obj$summary$sample_var_composition["Polygenic"] <- var(polygenic)
  obj$summary$expected_total_var <- sum(obj$summary$expected_var_composition)
  obj$summary$sample_total_var <- var(obj$data$y)

  # Calculate covariance v
  obj$covariance$V <- covariance(obj, calculate = TRUE)

  obj
}
