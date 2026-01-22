#' Reconstruct scenario-level multipliers used to build W
#'
#' This matches the logic in est_base_weight_mat(): we expand
#' ratios, and (optionally) ratiostau and ratiosalpha into a list
#' of S scenarios, where S = dim(W)[3].
#'
#' @export
build_scenario_table <- function(object) {
  ratios      <- object@ratios
  ratiostau   <- object@ratiostau
  ratiosalpha <- object@ratiosalpha
  subgroup    <- object@subgroup
  tvp         <- object@tvp
  
  out <- list()
  
  if (!subgroup && !tvp) {
    # simplest case: one dimension = ratios
    for (r in seq_along(ratios)) {
      out[[length(out) + 1]] <- data.frame(
        scenario    = length(out) + 0,
        ratio       = ratios[r],
        ratiostau   = NA_real_,
        ratiosalpha = NA_real_
      )
    }
  } else if (subgroup && !tvp) {
    for (r in seq_along(ratios)) {
      for (t in seq_along(ratiostau)) {
        out[[length(out) + 1]] <- data.frame(
          scenario    = length(out) + 0,
          ratio       = ratios[r],
          ratiostau   = ratiostau[t],
          ratiosalpha = NA_real_
        )
      }
    }
  } else if (tvp && !subgroup) {
    for (r in seq_along(ratios)) {
      for (a in seq_along(ratiosalpha)) {
        out[[length(out) + 1]] <- data.frame(
          scenario    = length(out) + 0,
          ratio       = ratios[r],
          ratiostau   = NA_real_,
          ratiosalpha = ratiosalpha[a]
        )
      }
    }
  } else {
    # (subgroup == TRUE && tvp == TRUE) -- rare but possible
    for (r in seq_along(ratios)) {
      for (t in seq_along(ratiostau)) {
        for (a in seq_along(ratiosalpha)) {
          out[[length(out) + 1]] <- data.frame(
            scenario    = length(out) + 0,
            ratio       = ratios[r],
            ratiostau   = ratiostau[t],
            ratiosalpha = ratiosalpha[a]
          )
        }
      }
    }
  }
  
  do.call(rbind, out)
}
