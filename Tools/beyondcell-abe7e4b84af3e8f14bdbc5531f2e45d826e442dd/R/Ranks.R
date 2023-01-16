#' @title Ranks the signatures from most sensitive to least sensitive
#' @description  This function computes the beyondcell score's (BCS) statistics
#' of each signature and ranks them according to the switch point and mean.
#' @name bcRanks
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column of interest. If
#' \code{idents = NULL}, the function computes the ranks using all cells. If
#' \code{idents != NULL}, the signatures' ranks are computed for each level in
#' \code{idents}.
#' @param extended If \code{extended = TRUE}, this function returns the switch
#' point, mean, median, standard deviation, variance, min, max, proportion of
#' \code{NaN} and residuals' mean per signature. If \code{extended = FALSE},
#' this function returns only the switch point, mean and residuals' mean.
#' @return A \code{beyondcell} object with the results in a new entry of
#' \code{bc@@ranks}: \code{bc@@ranks[["general"]]} (if \code{idents = NULL}) or
#' \code{bc@@ranks[[idents]]} (if \code{idents != NULL}).
#' @examples
#' @export

bcRanks <- function(bc, idents = NULL, extended = TRUE) {
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop('Idents must be a single metadata column.')
    }
    if (idents %in% colnames(bc@meta.data)) {
      if (idents %in% names(bc@ranks)) {
        warning(paste0('$', idents, ' already exists in bc@ranks. ',
                       'Entry has been overwritten.'))
      }
      meta <- bc@meta.data[colnames(bc@normalized), idents, drop = TRUE]
    } else {
      stop('Idents not found.')
    }
  } else {
    if ("general" %in% names(bc@ranks)) {
      warning('$general already exists in bc@ranks. Entry has been overwritten.')
    }
  }
  # Check extended.
  if (length(extended) != 1 | !is.logical(extended[1])) {
    stop('extended must be TRUE or FALSE.')
  }
  # --- Code ---
  # Signatures in bc.
  sigs <- rownames(bc@normalized)
  n.rows <- length(sigs)
  # Cells in bc.
  cells <- colnames(bc@normalized)
  # Extended.
  if (extended) n <- 6
  else n <- 1
  # Statistics for all normalized BCS.
  if (is.null(idents)) {
    # Progress bar.
    total <- n.rows * n
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    # Name of the output.
    idents <- "general"
    # Column to order by.
    order.col <- "rank"
    # Get Statistics dataframe.
    final.stats <- GetStatistics(bc = bc, signatures = sigs, cells = cells,
                                 pb = pb, total = total, i = 1, n.rows = n.rows,
                                 extended = extended)
    # Add rank.
    sig.order <- order(-1 * final.stats$switch.point, final.stats$mean,
                       decreasing = TRUE)
    final.stats$rank[sig.order] <- 1:nrow(final.stats)
    # Reorder columns.
    final.stats <- final.stats[c("rank", colnames(final.stats)[-ncol(final.stats)])]
  } else {
    # Metadata levels.
    lvls <- sort(unique(as.factor(meta)), decreasing = FALSE)
    # Progress bar.
    total <- n.rows * n * length(lvls)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    # Column to order by.
    order.col <- paste0("rank.", lvls[1])
    # For each level...
    stats <- lapply(seq_along(lvls), function(i) {
      ### Cells belonging to that level.
      group.cells <- rownames(bc@meta.data)[row.names(bc@meta.data) %in% cells &
                                              bc@meta.data[, idents] == lvls[i]]
      group.cells <- group.cells[!is.na( group.cells)]
      ### Subset bc with group.cells.
      sub.bc <- bc
      sub.bc@regression <- list(order = rep("", 2), vars = NULL,
                                order.background = rep("", 2))
      sub.bc@background <- matrix(ncol = 0, nrow = 0)
      sub.bc <- suppressWarnings(suppressMessages(
        bcSubset(sub.bc, cells = group.cells)))
      ### Get Statistics dataframe.
      out <- GetStatistics(bc = sub.bc, signatures = sigs, cells = group.cells,
                           pb = pb, total = total, i = i, n.rows = n.rows,
                           extended = extended)
      ### Add 4squares group annotation.
      # Get info about drugs (their corresponding name in bc, the preferred name
      # used by beyondcell and the MoA).
      info <- FindDrugs(sub.bc, x = rownames(sub.bc@scaled), na.rm = FALSE)
      # Get dataset switch points.
      sp <- data.frame(switch.point = bc@switch.point[info$bc.names],
                       row.names = info$bc.names)
      # Subset residuals' means and switch points.
      res <- out[, "residuals.mean", drop = FALSE]
      df <- transform(merge(res, sp, by = 0), row.names = Row.names, Row.names = NULL)
      # Residual's deciles.
      res.decil <- quantile(as.numeric(out$residuals.mean), 
                            prob = seq(from = 0, to = 1, length = 11))
      # Group annotation.
      sp_lower_01 <- as.numeric(df$switch.point) < 0.1
      sp_lower_06 <- as.numeric(df$switch.point) < 0.6
      sp_higher_04 <- as.numeric(df$switch.point) > 0.4
      sp_higher_09 <- as.numeric(df$switch.point) > 0.9
      res_lower_10 <- as.numeric(df$residuals.mean) < res.decil[["10%"]]
      res_higher_90 <- as.numeric(df$residuals.mean) > res.decil[["90%"]]
      df$group <- rep("", times = nrow(df))
      df$group[sp_lower_01 & res_higher_90] <- "TOP-HighSensitivity"
      df$group[sp_higher_09 & res_lower_10] <- "TOP-LowSensitivity"
      df$group[sp_higher_04 & sp_lower_06 & res_lower_10] <- "TOP-Differential-LowSensitivity"
      df$group[sp_higher_04 & sp_lower_06 & res_higher_90] <- "TOP-Differential-HighSensitivity"
      # Add to dataframe.
      out$group <- df$group
      ### Add rank.
      sig.order <- order(-1 * out$switch.point, out$mean, decreasing = TRUE)
      out$rank[sig.order] <- 1:nrow(out)
      ### Reorder columns.
      out <- out[c("rank", colnames(out)[-ncol(out)])]
      ### Add level suffix to colnames.
      colnames(out) <- paste0(colnames(out), ".", lvls[i])
      return(out)
    })
    # Merge dataframes of all levels.
    final.stats <- do.call(what = cbind.data.frame, args = stats)
  }
  # Add Drug name and MoA to final.stats.
  cols <- colnames(final.stats)
  info <- subset(drugInfo, subset = drugInfo$IDs %in% rownames(final.stats))
  if (dim(info)[1] > 0) {
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(x) {
      paste(na.omit(unique(x)), collapse = "; ")
    })
  }
  rownames(info) <- info$IDs
  info <- info[, c("drugs", "preferred.drug.names", "MoAs", "targets", "sources")]
  final.stats <- transform(merge(final.stats, info, by = 0, all.x = TRUE),
                           row.names = Row.names, Row.names = NULL)
  # Order by rank and reorder columns.
  final.stats <- final.stats[order(final.stats[, order.col], decreasing = FALSE),
                             c("drugs", "preferred.drug.names", "MoAs",
                               "targets", "sources", cols)]
  # Add to beyondcell object.
  bc@ranks[[idents]] <- final.stats
  # Close the progress bar.
  Sys.sleep(0.1)
  close(pb)
  return(bc)
}

#' @title Returns the first/last n ranked signatures
#' @description  This function returns the top/bottom \code{n} signatures ranked
#' by \code{\link[beyondcell]{bcRanks}}. If the rank has not been previously
#' computed, \code{rankSigs} performs the ranking itself.
#' @name rankSigs
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column of interest. If
#' \code{idents = NULL}, the function uses the general rank computed with all
#' cells.
#' @param cond Level of \code{idents} to rank by the output vector. If
#' \code{idents = NULL}, this parameter is deprecated.
#' @param n Number of signatures to return in the output vector.
#' @param decreasing Logical. Return the top \code{n} signatures (default) or
#' the bottom \code{n} signatures (\code{decreasing = FALSE})?.
#' @return An ordered vector with the signature's names.
#' @examples
#' @export

rankSigs <- function(bc, idents = NULL, cond = NULL, n = 10,
                     decreasing = TRUE) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop('Idents must be a single metadata column.')
    }
    if (!idents %in% colnames(bc@meta.data)) {
      stop('Idents not found.')
    }
    if (is.null(cond[1])) {
      stop('Invalid cond.')
    }
    meta <- idents
  } else {
    meta <- "general"
    if (!is.null(cond[1])) {
      warning('idents not specified, cond is deprecated.')
      cond <- NULL
    }
  }
  # Check cond.
  if (!is.null(cond)) {
    if (length(cond) != 1) {
      stop('cond must be a single idents level')
    }
    if (!cond %in% levels(as.factor(bc@meta.data[, idents]))) {
      stop(paste0(cond, ' is not a level of ', idents, '.'))
    }
  }
  # Check n.
  if (length(n) != 1 | (!is.numeric(n) & !is.character(n))) {
    stop('n must be a single number or "all".')
  }
  if (is.numeric(n) & (n < 1 | n%%1 != 0)) stop('n must be an integer > 0.')
  else if (is.character(n)) {
    if (n == "all") n <- nrow(bc@normalized)
    else stop('To select all signatures, please set n = "all".')
  }
  # Check decreasing.
  if (length(decreasing) != 1 | !is.logical(decreasing)) {
    stop('decreasimg must be TRUE or FALSE.')
  }
  # --- Code ---
  # If ranks have not been computed, compute them now.
  if (!meta %in% names(bc@ranks)) {
    message('Computing ranks...')
    bc <- bcRanks(bc, idents = idents, extended = FALSE)
  }
  # Get ranks for the specified idents.
  df <- bc@ranks[[meta]]
  # Ranks to select.
  if (decreasing == TRUE) {
    idx <- 1:n
  } else {
    idx <- nrow(bc@normalized):(nrow(bc@normalized) - n + 1)
  }
  # Return signatures whose rank == idx.
  order.col <- ifelse(test = meta == "general", yes = "rank",
                      no = paste0("rank.", cond))
  ordered.df <- df[order(df[, order.col], decreasing = FALSE), ]
  sigs <- ordered.df[idx, "Name"]
  return(sigs)
}
