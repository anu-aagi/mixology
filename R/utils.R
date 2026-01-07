#' @import Matrix
NULL

#' @import dplyr
NULL

.check.vcov.matrix <- function(size, vcov) {
  m <- NULL

  if (!(inherits(vcov, "Matrix") || inherits(vcov, "matrix"))) {
    m <- "Not a matrix or a Matrix object."
  } else if (!all(dim(vcov) == size)) {
    m <- "Incorrect dimension."
  } else if (!(is.numeric(vcov) || inherits(vcov, "Matrix")) || anyNA(vcov)) {
    m <- "Not numeric or contain NA values."
  } else if (!all.equal(vcov, t(vcov), tolerance = 1e-10)) {
    m <- "Not symmetric."
  }

  return(m)
}

.capture_all_args <- function(fun, env = parent.frame(), drop = NULL) {
  fmls <- formals(fun)

  args <- lapply(names(fmls), function(nm) {
    get(nm, envir = env)
  })
  names(args) <- names(fmls)

  args <- args[setdiff(names(args), drop)]
  args
}
