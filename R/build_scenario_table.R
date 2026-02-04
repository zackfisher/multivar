#' Reconstruct scenario-level multipliers used to build W
#'
#' This matches the logic in est_base_weight_mat(): we expand
#' ratios_unique, and (optionally) ratios_subgroup and ratios_unique_tvp into a list
#' of S scenarios, where S = dim(W)[3].
#'
#' @export
build_scenario_table <- function(object) {
  ratios_unique      <- object@ratios_unique
  ratios_subgroup   <- object@ratios_subgroup
  ratios_unique_tvp <- object@ratios_unique_tvp
  subgroup    <- object@subgroup
  tvp         <- object@tvp
  
  out <- list()
  
  if (!subgroup && !tvp) {
    # simplest case: one dimension = ratios_unique
    for (r in seq_along(ratios_unique)) {
      out[[length(out) + 1]] <- data.frame(
        scenario    = length(out) + 0,
        ratios_unique       = ratios_unique[r],
        ratios_subgroup   = NA_real_,
        ratios_unique_tvp = NA_real_
      )
    }
  } else if (subgroup && !tvp) {
    for (r in seq_along(ratios_unique)) {
      for (t in seq_along(ratios_subgroup)) {
        out[[length(out) + 1]] <- data.frame(
          scenario    = length(out) + 0,
          ratios_unique       = ratios_unique[r],
          ratios_subgroup   = ratios_subgroup[t],
          ratios_unique_tvp = NA_real_
        )
      }
    }
  } else if (tvp && !subgroup) {
    for (r in seq_along(ratios_unique)) {
      for (a in seq_along(ratios_unique_tvp)) {
        out[[length(out) + 1]] <- data.frame(
          scenario    = length(out) + 0,
          ratios_unique       = ratios_unique[r],
          ratios_subgroup   = NA_real_,
          ratios_unique_tvp = ratios_unique_tvp[a]
        )
      }
    }
  } else {
    # (subgroup == TRUE && tvp == TRUE) -- rare but possible
    for (r in seq_along(ratios_unique)) {
      for (t in seq_along(ratios_subgroup)) {
        for (a in seq_along(ratios_unique_tvp)) {
          out[[length(out) + 1]] <- data.frame(
            scenario    = length(out) + 0,
            ratios_unique       = ratios_unique[r],
            ratios_subgroup   = ratios_subgroup[t],
            ratios_unique_tvp = ratios_unique_tvp[a]
          )
        }
      }
    }
  }
  
  do.call(rbind, out)
}
