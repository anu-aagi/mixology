#' Simulate SNP data
#' @export
snp_simulation <- function(N = 300,
                           M = 600,
                           K = 10,
                           J = 1,
                           p_f = NULL, # population probabilities
                           bias = 1, # alpha scaling (shape)
                           variance = 1, # alpha scaling (spread)
                           rho = 0.99, # AR(1) autocorrelation
                           p_AA = 0.5, # DTMC: P(1 -> 1)
                           p_CC = p_AA, # DTMC: P(-1 -> -1)
                           sigma_l = 2,
                           sigma_lk = NULL, # Var(Amn) ~ sigma2_l + 1/K sigma2_lk + sigma2_e ~ 3.1
                           sigma_e = 1,
                           p_linkage = 0.99,
                           p_mate = 1 / J,
                           seed = NULL) {
  # Get call
  args <- .capture_all_args(snp_simulation)

  if (is.numeric(seed)) set.seed(seed)

  # Check population probabilities p_f
  if (!is.null(p_f)) {
    if (!is.numeric(p_f) || any(p_f < 0)) {
      stop("Bad p_f specification")
    } else if (sum(p_f) != 1) {
      message("Scale p_f to sum to 1")
      p_f <- p_f / sum(p_f)
    }
    message("Suppress J by length(p_f)")
    J <- length(p_f)
  } else {
    p_f <- rep(1, J) / J
  }

  # Check sigma_lk and possibly override K
  if (is.null(sigma_lk)) {
    p <- (seq_len(K) / (K + 1))^bias
    p <- p / sum(p)
    sigma_lk <- sqrt(rep(1 / sum(p^2), K))
  } else if (!is.numeric(sigma_lk) || any(sigma_lk < 0)) {
    stop("Bad sigma_lk specification")
  } else if (length(sigma_lk) == 1) {
    sigma_lk <- rep(sigma_lk, K)
  } else {
    message("Suppress K by length(sigma_lk)")
    K <- length(sigma_lk)
  }


  # Simulate factor matrix f_mat and population labels z_f
  variance <- variance
  bias <- bias

  alpha <- replicate(
    K,
    (1 / variance) * runif(J)^bias
  )
  alpha <- matrix(alpha, nrow = J) # J-by-K
  alpha[alpha > 1e16] <- 1e16
  alpha[alpha < 1e-16] <- 1e-16

  f_mat <- matrix(0, nrow = N, ncol = K)
  z_f <- sample(seq_len(J), N, replace = TRUE, prob = p_f)

  # Check if we are going to simulating mating event
  if (!is.null(p_mate)) {
    do_mate <- TRUE
    f_mat_partner <- f_mat
    z_f_partner <- numeric(N)
  } else {
    do_mate <- FALSE
    f_mat_partner <- z_f_partner <- NULL
  }

  # Self mixing
  for (j in seq_len(J)) {
    idx <- z_f == j
    alpha_j <- alpha[j, ]
    n.j <- sum(idx)
    if (any(idx)) {
      f_mat[idx, ] <- MCMCpack::rdirichlet(n.j, alpha_j)
      if (do_mate) {
        z_f_partner[idx] <- sample(
          c(j, seq_len(J)[-j]),
          n.j,
          replace = TRUE,
          prob = c(p_mate, rep((1 - p_mate) / (J - 1), J - 1))
        )
      }
    }
  }

  # Partner mixing
  if (do_mate) {
    for (j in seq_len(J)) {
      idx <- z_f_partner == j
      alpha_j <- alpha[j, ]
      n.j <- sum(idx)
      if (any(idx)) {
        f_mat_partner[idx, ] <- MCMCpack::rdirichlet(n.j, alpha_j)
      }
    }
  }

  # Name factor-related objects
  z_f <- paste0("Population", z_f) |>
    factor(levels = paste0("Population", 1:J))
  rownames(f_mat) <- names(z_f) <- paste0("Sample", 1:N)
  colnames(f_mat) <- colnames(alpha) <- paste0("Factor", 1:K)
  rownames(alpha) <- levels(z_f)
  if (do_mate) {
    dimnames(f_mat_partner) <- dimnames(f_mat)
    z_f_partner <- paste0("Population", z_f_partner) |>
      factor(levels = levels(z_f))
  }

  # Helper: simulate one AR(1) loading profile
  sim_one_ar1 <- function(sigma, rho, p_linkage, M) {
    u <- numeric(M)
    z <- numeric(M)
    u[1] <- rnorm(1, 0, sigma)
    for (m in 2:M) {
      z <- sample(c(1, 0), 1, prob = c(p_linkage, 1 - p_linkage))
      u[m] <- z * (rho * u[m - 1] + rnorm(1, 0, sigma * sqrt(1 - rho^2)))
    }
    u
  }

  # Helper: simulate two-state DTMC (+1 / -1)
  sim_one_dtmc <- function(p_AA, p_CC, M) {
    z <- numeric(M)
    p_ini <- (1 - p_CC) / (2 - p_AA - p_CC)
    z[1] <- sample(c(1, -1), 1, prob = c(p_ini, 1 - p_ini))
    for (m in 2:M) {
      if (z[m - 1] == 1) {
        z[m] <- sample(c(1, -1), 1, prob = c(p_AA, 1 - p_AA))
      } else {
        z[m] <- sample(c(1, -1), 1, prob = c(1 - p_CC, p_CC))
      }
    }
    z
  }

  # Simulate shared and loading-specific AR(1) processes
  u_lk <- sapply(sigma_lk, sim_one_ar1, rho = rho, p_linkage = p_linkage, M = M)
  u_lk <- matrix(u_lk, nrow = M)
  u_l <- sim_one_ar1(sigma_l, rho, p_linkage, M)

  # Simulate shared and loading-specific DTMCs for allele preference
  z_lk <- replicate(K, sim_one_dtmc(p_AA, p_CC, M), simplify = "vector")
  z_lk <- matrix(z_lk, nrow = M)
  z_l <- sim_one_dtmc(p_AA, p_CC, M)

  # Combine into loading matrix l_mat
  l_mat <- abs(u_lk) * z_lk + abs(u_l) * z_l

  # Name marker-related objects
  rownames(u_lk) <- rownames(z_lk) <-
    names(u_l) <- names(z_l) <-
    rownames(l_mat) <- paste0("Marker", 1:M)
  colnames(u_lk) <- colnames(z_lk) <- colnames(l_mat) <- paste0("Factor", 1:K)

  # Genotype intensity matrix A
  A <- tcrossprod(f_mat, l_mat) +
    matrix(rnorm(N * M, sd = sigma_e), nrow = N)
  A_partner <- if (do_mate) {
    tcrossprod(f_mat_partner, l_mat) +
      matrix(rnorm(N * M, sd = sigma_e), nrow = N)
  } else {
    NULL
  }

  # Binomial sampling of genotypes (sum of two Bernoulli draws)
  X_prob <- pnorm(A)
  X <- matrix(runif(N * M), nrow = N) < X_prob
  if (do_mate) X_prob <- pnorm(A_partner)
  X <- X + (matrix(runif(N * M), nrow = N) < X_prob)
  X <- Matrix::Matrix(X)

  # Return simulation components
  obj <- list(
    data = list(
      y       = NULL,
      X       = X,
      W       = NULL,
      N       = N,
      M       = M,
      J       = J,
      meta    = list(
        Membership =  data.frame(Membership1 = z_f, Membership2 = z_f_partner),
        Environment = NULL,
        Misc = NULL
      ),
      added = NULL,
      samples = rownames(X),
      markers = colnames(X)
    ),
    design = list(
      fixed = list(Additive = NULL),
      random = list(
        Polygenic = NULL,
        Environment = NULL
      )
    ),
    effects = list(
      fixed = list(Additive = NULL),
      random = list(
        Polygenic = NULL,
        Environment = NULL
      )
    ),
    covariance = list(
      V = NULL,
      R = NULL,
      G = list(
        Polygenic = NULL,
        Environment = NULL
      )
    ),
    component = NULL,
    summary = list(
      mean = NULL,
      expected_total_var = NULL,
      sample_total_var = NULL,
      expected_var_composition = NULL,
      sample_var_composition = NULL
    ),
    params = list(
      snp_simulation = args
    ),
    misc = list(
      snp_simulation = list(
        A             = A,
        A_partner     = A_partner,
        alpha         = alpha,
        f_mat         = f_mat,
        f_mat_partner = f_mat_partner,
        l_mat         = l_mat,
        u_lk          = u_lk,
        z_lk          = z_lk,
        u_l           = u_l,
        z_l           = z_l
      )
    ),
    seed = seed
  )

  # Return the object
  class(obj) <- c("cocktail", "list")

  obj
}
