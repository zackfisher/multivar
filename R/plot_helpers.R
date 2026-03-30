# ==============================================================================
# Multivar Plotting Helper Functions
# ==============================================================================
#
# Internal functions for plotting transition matrices.
# These are not exported and serve as building blocks for user-facing plots.
#

#' Core transition matrix heatmap plotting engine
#'
#' @param mat_list List of matrices to plot, or a single matrix
#' @param titles Character vector of titles for each matrix (optional)
#' @param facet_ncol Number of columns for faceting (when plotting multiple)
#' @param lb Lower bound for color scale
#' @param ub Upper bound for color scale
#' @param show_zeros Logical. If FALSE (default), zeros are shown as white/NA
#' @param palette Character. Color palette to use ("default", "viridis", "greyscale")
#' @param show_dimnames Logical. Show row/column names?
#' @param title Main title for the plot
#' @param subtitle Subtitle for the plot
#'
#' @return A ggplot2 object
#' @keywords internal
.plot_transition_heatmap <- function(mat_list,
                                     titles = NULL,
                                     facet_ncol = 3,
                                     lb = -1,
                                     ub = 1,
                                     show_zeros = FALSE,
                                     palette = "default",
                                     show_dimnames = TRUE,
                                     title = NULL,
                                     subtitle = NULL) {

  # Check dependencies
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it with: install.packages('ggplot2')",
         call. = FALSE)
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is required for plotting. Please install it with: install.packages('reshape2')",
         call. = FALSE)
  }

  # Convert single matrix to list
  if (is.matrix(mat_list)) {
    mat_list <- list(mat_list)
  }

  # Validate input
  if (!is.list(mat_list)) {
    stop("mat_list must be a matrix or list of matrices")
  }

  # Generate default titles if needed
  if (is.null(titles)) {
    if (length(mat_list) == 1) {
      titles <- ""
    } else {
      titles <- paste0("Matrix ", seq_along(mat_list))
    }
  }

  # Ensure variable names exist
  for (i in seq_along(mat_list)) {
    mat <- mat_list[[i]]
    if (is.null(colnames(mat))) {
      colnames(mat) <- paste0("V", seq_len(ncol(mat)))
    }
    if (is.null(rownames(mat))) {
      rownames(mat) <- paste0("V", seq_len(nrow(mat)))
    }
    mat_list[[i]] <- mat
  }

  # Convert matrices to long-form dataframe
  df_list <- lapply(seq_along(mat_list), function(i) {
    mat <- mat_list[[i]]
    df <- stats::setNames(reshape2::melt(mat), c('rows', 'vars', 'values'))
    df$panel <- titles[i]
    df
  })

  df <- do.call("rbind", df_list)

  # Handle zeros
  if (!show_zeros) {
    df$values[df$values == 0] <- NA
  }

  # Set up color palette
  if (palette == "default") {
    zf_red <- grDevices::rgb(255, 0, 90, maxColorValue = 255)
    zf_blue <- grDevices::rgb(0, 152, 233, maxColorValue = 255)
    zf_fore <- "white"

    colfunc_low <- grDevices::colorRampPalette(c(zf_red, zf_fore))
    colfunc_high <- grDevices::colorRampPalette(c(zf_fore, zf_blue))
    colors_to_use <- c(colfunc_low(6)[1:3], zf_fore, colfunc_high(6)[4:6])

  } else if (palette == "viridis") {
    if (!requireNamespace("viridis", quietly = TRUE)) {
      warning("Package 'viridis' not available. Using default palette.")
      palette <- "default"
      # Recursively call with default
      return(.plot_transition_heatmap(mat_list, titles, facet_ncol, lb, ub,
                                      show_zeros, palette = "default",
                                      show_dimnames, title, subtitle))
    }
    colors_to_use <- NULL  # Will use scale_fill_viridis instead

  } else if (palette == "greyscale") {
    colors_to_use <- c("black", "grey30", "grey60", "white", "grey60", "grey30", "black")

  } else {
    stop("Unknown palette: ", palette, ". Use 'default', 'viridis', or 'greyscale'")
  }

  # Theme colors
  text_color <- grDevices::rgb(51, 51, 51, maxColorValue = 255)
  plot_background <- "white"

  # Set factor levels for correct ordering
  var_levels <- colnames(mat_list[[1]])
  row_levels <- rev(rownames(mat_list[[1]]))

  df$rows <- factor(df$rows, levels = row_levels)
  df$vars <- factor(df$vars, levels = var_levels)

  # Color scale limits
  limit <- max(abs(c(lb, ub))) * c(-1, 1)

  # Build ggplot
  gg <- ggplot2::ggplot(df, ggplot2::aes(y = .data$rows, 
                                         x = .data$vars, 
                                         fill = .data$values))
  gg <- gg + ggplot2::geom_tile()

  # Apply color scale
  if (palette == "viridis") {
    gg <- gg + viridis::scale_fill_viridis(
      option = "magma",
      limits = limit,
      na.value = "white",
      guide = ggplot2::guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        ticks.linewidth = 1,
        frame.linewidth = 1
      )
    )
  } else {
    gg <- gg + ggplot2::scale_fill_gradientn(
      colors = colors_to_use,
      limits = limit,
      na.value = "white",
      guide = ggplot2::guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        ticks.linewidth = 1,
        frame.linewidth = 1
      )
    )
  }

  # Apply theme
  gg <- gg + ggplot2::coord_equal()
  gg <- gg + ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(hjust = 0, color = text_color, face = "bold"),
    strip.text = ggplot2::element_text(hjust = 0, color = text_color, size = 12, face = "bold"),
    strip.background = ggplot2::element_rect(fill = plot_background, color = plot_background),
    panel.spacing.x = grid::unit(0.5, "cm"),
    panel.spacing.y = grid::unit(0.5, "cm"),
    legend.background = ggplot2::element_rect(fill = plot_background, color = plot_background),
    legend.title = ggplot2::element_text(size = 12, color = text_color),
    legend.title.align = 1,
    legend.text = ggplot2::element_text(size = 10, color = text_color),
    legend.text.align = 1,
    plot.background = ggplot2::element_rect(fill = plot_background, color = plot_background),
    panel.border = ggplot2::element_rect(fill = NA, colour = 'black', linewidth = 1)
  )

  # Handle axis text
  if (show_dimnames) {
    gg <- gg + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, color = text_color),
      axis.text.y = ggplot2::element_text(size = 10, color = text_color)
    )
  } else {
    gg <- gg + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  }

  # Labels
  gg <- gg + ggplot2::labs(fill = '', x = NULL, y = NULL)

  if (!is.null(title)) {
    gg <- gg + ggplot2::labs(title = title)
  }

  if (!is.null(subtitle)) {
    gg <- gg + ggplot2::labs(subtitle = subtitle)
  }

  # Faceting for multiple panels
  if (length(mat_list) > 1) {
    gg <- gg + ggplot2::facet_wrap(~ panel, ncol = facet_ncol)
  }

  return(gg)
}


#' Extract matrices from a multivar fit object
#'
#' @param fit Result from cv.multivar() or similar
#' @param type Type of matrices to extract ("common", "unique", "total", "subgrp", "intercept", "tvp", "common_tvp")
#' @param subjects Which subjects to include ("all" or numeric vector)
#' @param periods Which periods to include for TVP models ("all" or numeric vector)
#'
#' @return List of matrices with appropriate names
#' @keywords internal
.extract_matrices_from_fit <- function(fit,
                                       type = c("common", "unique", "total", "subgrp", "intercept", "tvp", "common_tvp"),
                                       subjects = "all",
                                       periods = "all") {

  type <- match.arg(type)

  # Access the matrices
  mats <- fit$mats

  if (type == "common") {
    if (is.null(mats$common)) {
      stop("No common effects in this model (K=1 or model without common effects)")
    }
    return(list("Common Effects" = mats$common))

  } else if (type == "unique") {
    if (is.null(mats$unique)) {
      stop("No unique effects in this model")
    }

    # Select subjects
    if (subjects[1] == "all") {
      subject_idx <- seq_along(mats$unique)
    } else {
      subject_idx <- subjects
    }

    mat_list <- mats$unique[subject_idx]
    names(mat_list) <- paste0("Subject ", subject_idx)
    return(mat_list)

  } else if (type == "total") {
    if (is.null(mats$total)) {
      stop("No total effects in this model")
    }

    # Select subjects
    if (subjects[1] == "all") {
      subject_idx <- seq_along(mats$total)
    } else {
      subject_idx <- subjects
    }

    mat_list <- mats$total[subject_idx]
    names(mat_list) <- paste0("Subject ", subject_idx)
    return(mat_list)

  } else if (type == "subgrp") {
    if (is.null(mats$subgrp)) {
      stop("No subgroup effects in this model (K=1 or model without subgroups)")
    }

    # Subgroup effects - return all unique subgroup matrices
    # Note: subgrp may be a list (one per subject) or list of unique subgroup matrices
    if (is.list(mats$subgrp)) {
      # Get unique subgroup matrices
      unique_subgrps <- unique(mats$subgrp)
      mat_list <- unique_subgrps
      names(mat_list) <- paste0("Subgroup ", seq_along(unique_subgrps))
    } else {
      mat_list <- list("Subgroup Effects" = mats$subgrp)
    }
    return(mat_list)

  } else if (type == "intercept") {
    if (is.null(mats$intercepts)) {
      stop("No intercepts in this model (fit without intercept=TRUE)")
    }

    # Select subjects
    if (subjects[1] == "all") {
      subject_idx <- seq_along(mats$intercepts)
    } else {
      subject_idx <- subjects
    }

    # Convert intercept vectors to matrices for plotting
    mat_list <- lapply(subject_idx, function(i) {
      as.matrix(mats$intercepts[[i]])
    })
    names(mat_list) <- paste0("Subject ", subject_idx, " Intercept")
    return(mat_list)

  } else if (type == "tvp") {
    if (is.null(mats$tvp)) {
      stop("No TVP effects in this model (fit without tvp=TRUE)")
    }

    # Select subjects
    if (subjects[1] == "all") {
      subject_idx <- seq_along(mats$tvp)
    } else {
      subject_idx <- subjects
    }

    # TVP structure depends on K
    # For K=1: mats$tvp[[1]] is a list of period matrices
    # For K>1: mats$tvp[[i]] is a list of time-point matrices

    k <- length(mats$tvp)

    if (k == 1) {
      # K=1: TVP effects are organized by period
      period_mats <- mats$tvp[[1]]

      # Select periods
      if (periods[1] == "all") {
        period_idx <- seq_along(period_mats)
      } else {
        period_idx <- periods
      }

      mat_list <- period_mats[period_idx]
      names(mat_list) <- paste0("TVP Period ", period_idx)
      return(mat_list)

    } else {
      # K>1: TVP effects organized by subject, then by time point
      # For visualization, we'll extract by periods (assuming breaks define periods)

      # Check if we have period structure
      if (is.list(mats$tvp[[subject_idx[1]]])) {
        # Get number of periods from first subject
        num_periods <- length(mats$tvp[[subject_idx[1]]])

        # Select periods
        if (periods[1] == "all") {
          period_idx <- seq_len(num_periods)
        } else {
          period_idx <- periods
        }

        # Extract matrices for selected subjects and periods
        mat_list <- list()
        for (s in subject_idx) {
          for (p in period_idx) {
            mat_list[[length(mat_list) + 1]] <- mats$tvp[[s]][[p]]
            names(mat_list)[length(mat_list)] <- paste0("Subject ", s, " Period ", p)
          }
        }
        return(mat_list)
      } else {
        stop("TVP structure not recognized. Expected nested list.")
      }
    }

  } else if (type == "common_tvp") {
    if (is.null(mats$common_tvp)) {
      stop("No common TVP effects in this model (fit without tvp=TRUE and common_tvp_effects=TRUE)")
    }

    # common_tvp is a list of period matrices
    # Select periods
    if (periods[1] == "all") {
      period_idx <- seq_along(mats$common_tvp)
    } else {
      period_idx <- periods
    }

    mat_list <- mats$common_tvp[period_idx]
    names(mat_list) <- paste0("Common TVP Period ", period_idx)
    return(mat_list)
  }
}


#' Extract matrices from a multivar simulation object
#'
#' @param sim Result from multivar_sim() or multivar_sim_subgroups()
#' @param type Type of matrices to extract ("common", "unique", "total", "subgrp")
#' @param subjects Which subjects to include ("all" or numeric vector)
#'
#' @return List of matrices with appropriate names
#' @keywords internal
.extract_matrices_from_sim <- function(sim,
                                       type = c("common", "unique", "total", "subgrp"),
                                       subjects = "all") {

  type <- match.arg(type)

  if (type == "common") {
    if (is.null(sim$mat_com)) {
      stop("No common effects in this simulation")
    }
    return(list("True Common Effects" = sim$mat_com))

  } else if (type == "unique") {
    if (is.null(sim$mat_ind_unique)) {
      stop("No unique effects in this simulation")
    }

    # Select subjects
    if (subjects[1] == "all") {
      subject_idx <- seq_along(sim$mat_ind_unique)
    } else {
      subject_idx <- subjects
    }

    mat_list <- sim$mat_ind_unique[subject_idx]
    names(mat_list) <- paste0("Subject ", subject_idx, " (True)")
    return(mat_list)

  } else if (type == "total") {
    if (is.null(sim$mat_ind_final)) {
      stop("No total effects in this simulation")
    }

    # Select subjects
    if (subjects[1] == "all") {
      subject_idx <- seq_along(sim$mat_ind_final)
    } else {
      subject_idx <- subjects
    }

    mat_list <- sim$mat_ind_final[subject_idx]
    names(mat_list) <- paste0("Subject ", subject_idx, " (True)")
    return(mat_list)

  } else if (type == "subgrp") {
    if (is.null(sim$mat_sub_unique)) {
      stop("No subgroup effects in this simulation (use multivar_sim_subgroups)")
    }

    # Get unique subgroup matrices
    unique_subgrps <- unique(sim$mat_sub_unique)
    mat_list <- unique_subgrps
    names(mat_list) <- paste0("Subgroup ", seq_along(unique_subgrps), " (True)")
    return(mat_list)
  }
}
