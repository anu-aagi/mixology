
#' @export
#' @method print cocktail
print.cocktail <- function(obj) {
  message("Enjoy your cocktail :) \n")

  if (!is.null(obj$data$X)) {
    message("Genotype data exists: ")
    message(
      "  ", N, " samples, ",
      M, " markers, ",
      J, " populations ",
      ifelse(is.null(obj$params$snp_simulation$p_mate), "without", "with"),
      " random mating."
    )
  } else {
    message("Genotype data is absent. ")
  }

  message("\n")

  if (!is.null(obj$data$y)) {
    message("Trait data exists: ")
    message(
      "  Trait Mean:\n",
      "    ", round(obj$summary$mean["Expected"], 2), " (Expected); ",
      round(obj$summary$mean["Sample"], 2), " (Sample)"
    )

    message("  Expected Variance Composition (%): ")
    v <- c(
      Total = obj$summary$expected_total_var,
      obj$summary$expected_var_composition
    )
    v_prop <- v / v["Total"]
    v_prop <- round(v_prop * 100, digits = 2)
    message(
      "    ",
      paste0(names(v_prop), " ",
             v, " ",
             "(", v_prop, "%)",
             collapse = "; "
      )
    )

    message("  Sample Variance Composition (%): ")
    v <- c(
      Total = obj$summary$sample_total_var,
      obj$summary$sample_var_composition
    )
    v_prop <- v / v["Total"]
    v_prop <- round(v_prop * 100, digits = 2)
    message(
      "    ",
      paste0(names(v_prop), " ",
             round(v, 2), " ",
             "(", v_prop, "%)",
             collapse = "; "
      )
    )
  } else {
    message("Trait data is absent.")
  }
}

#' @export
#' @method dim cocktail
dim.cocktail <- function(obj) {
  c(obj$data$N, obj$data$M)
}

#' @export
#' @method dimnames cocktail
dimnames.cocktail <- function(obj) {
  list(obj$data$samples, obj$data$markers)
}
