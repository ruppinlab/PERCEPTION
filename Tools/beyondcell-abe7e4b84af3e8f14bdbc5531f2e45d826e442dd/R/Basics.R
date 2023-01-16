#' @title Substracts to the first element of a vector the rest of elements
#' @description This function substracts to the first element of a numeric
#' vector (\code{x[1]}) the rest of elements of the same vector.
#' (\code{x[2:length(x)]}).
#' @name minus
#' @param x Numeric vector.
#' @param na.rm (From \code{base}) Logical. Should missing values (including
#' \code{NaN}) be removed?
#' @return The result of the substraction.
#' @examples
#' @export

minus <- function(x, na.rm = FALSE) {
  # --- Checks ---
  # Check x.
  if (!is.numeric(x)) stop('x must be numeric.')
  # Check na.rm.
  if (length(na.rm) != 1 | !is.logical(na.rm)) {
    stop('na.rm must be TRUE or FALSE.')
  }
  # --- Code ---
  # If x is a single number, append a zero so we can run the next step.
  if (length(x) == 1) x <- c(x, 0)
  # Substract to the first element of x the rest of elements.
  out <- sum(x[1], na.rm = na.rm) - sum(x[2:length(x)], na.rm = na.rm)
  return(out)
}

#' @title Computes the column substraction
#' @description This function substracts to the first element of each column of
#' a rectangular object (\code{x[1, n]}) the rest of elements of the same column
#' (\code{x[2:length(x), n]}).
#' @name colMinus
#' @param x A matrix or a dataframe.
#' @param na.rm (From \code{base}) Logical. Should missing values (including
#' \code{NaN}) from rows \code{2:length(x)} be omitted from the calculations?
#' @return A numeric rectangular object with the result of the substraction.
#' @examples
#' @export

colMinus <- function(x, na.rm = FALSE) {
  # --- Checks ---
  # Check x.
  if (!is.matrix(x) & !is.data.frame(x)) {
    stop('x must be a matrix or a data.frame.')
  }
  # Check na.rm.
  if (length(na.rm) != 1 | !is.logical(na.rm)) {
    stop('na.rm must be TRUE or FALSE.')
  }
  # --- Code ---
  # If x has a single row, append a row of zeros so we can run the next step.
  if (dim(x)[1] == 1) x <- rbind(x, rep(0, times = ncol(x)))
  # Substract to the first row of x the rest of rows.
  first.row <- x[1, , drop = FALSE]
  out <- first.row - colSums(x[2:nrow(x), , drop = FALSE], na.rm = na.rm)
  return(out)
}

#' @title Computes the mean, the median and the sd of a vector
#' @description This function computes the mean, the median and the sd of a
#' vector.
#' @name Mean.Med.SD
#' @param x Numeric vector.
#' @param na.rm (From base) Logical. Should missing values (including NaN) be
#' removed?
#' @return A named numeric vector with the mean, median and sd of \code{x}.
#' @examples
#' @export

Mean.Med.SD <- function(x) {
  # --- Checks ---
  # Check that x is a numeric vector.
  if (length(x) <= 1 | !is.numeric(x)) {
    stop('x must be a numeric vector.')
  }
  # --- Code ---
  stats.mean <- mean(x, na.rm = TRUE)
  stats.median <- median(x, na.rm = TRUE)
  stats.sd <- sd(x, na.rm = TRUE)
  return(setNames(c(stats.mean, stats.median, stats.sd),
                  c("mean", "median", "sd")))
}

#' @title Computes the BCS' statistics and ranks
#' @description This function computes the beyondcell scores' (BCS) statistics
#' and ranks returned by \code{\link[beyondcell]{bcRanks}}.
#' @name GetStatistics
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param signatures Vector with the names of the signatures of interest.
#' @param cells Vector with the names of the cells of interest.
#' @param pb \code{\link[utils]{txtProgressBar}}.
#' @param total Number of iterations to complete the \code{pb}.
#' @param i Iteration number. Used to increase the \code{pb}.
#' @param n.rows Number of signatures. Used to increase the \code{pb}.
#' @param extended If \code{extended = TRUE}, this function returns the switch
#' point, mean, median, sd, variance, min, max, proportion of \code{NaN} and
#' residuals' mean per signature. If \code{extended = FALSE}, this function
#' returns only the switch point, mean and residuals' mean.
#'
#' @return A \code{data.frame} with the BCS' statistics and ranks.
#' @examples
#' @export

GetStatistics <- function(bc, signatures, cells, pb, total, i, n.rows,
                          extended) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check signatures.
  in.signatures <- !is.null(signatures) & signatures %in% rownames(bc@normalized)
  if (all(!in.signatures)) {
    stop('None of the specified signatures were found.')
  } else if (any(!in.signatures)) {
    warning(paste0('These signatures were not found in bc: ',
                   paste0(signatures[!in.signatures], collapse = ", "), '.'))
  }
  # Check cells.
  in.cells <- !is.null(cells) & cells %in% colnames(bc@normalized)
  if (all(!in.cells)) {
    stop('None of the specified cells were found.')
  } else if (any(!in.cells)) {
    warning(paste0('These cells were not found in bc: ',
                   paste0(cells[!in.cells], collapse = ", "), '.'))
  }
  # Check pb.
  if (class(pb) != "txtProgressBar") stop('pb must be a txtProgressBar object.')
  # Check total.
  if (length(total) != 1 | total[1] <= 0 | total[1]%%1 != 0) {
    stop('total must be a positive integer.')
  }
  # Check i.
  if (length(i) != 1 | i[1] <= 0 | i[1]%%1 != 0) {
    stop('i must be a positive integer.')
  }
  # Check n.rows.
  if (length(n.rows) != 1 | n.rows[1] <= 0 | n.rows[1]%%1 != 0) {
    stop('n.rows must be a positive integer.')
  }
  # Check extended.
  if (length(extended) != 1 | !is.logical(extended[1])) {
    stop('extended must be TRUE or FALSE.')
  }
  # --- Code ---
  # Signatures and cells.
  signatures <- signatures[in.signatures]
  cells <- cells[in.cells]
  # Bins: Number of iterations to complete each 1% of the progress bar.
  bins <- ceiling(total / 100)
  # Switch points per signature.
  switch.p <- bc@switch.point[signatures]
  if (extended) {
    # Dataframe with mean, median and sd per signature.
    data <- cbind(seq_len(n.rows) + (n.rows * 6) * (i - 1),
                  bc@data[signatures, cells])
    mean.med.sd <- as.data.frame(t(apply(data, 1, function(u) {
      mms <- round(Mean.Med.SD(u[-1]), digits = 2)
      ### Update the progress bar.
      if (u[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = u[1])
      }
      return(mms)
    })))
    # Variance per signature.
    data[, 1] <- data[, 1] + n.rows
    variance.bcscore <- apply(data, 1, function(v) {
      variance <- var(v[-1], na.rm = TRUE)
      ### Update the progress bar.
      if (v[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = v[1])
      }
      return(variance)
    })
    # Min normalized BCS per signature.
    data[, 1] <- data[, 1] + n.rows
    min.bcscore <- apply(data, 1, function(w) {
      min.bcs <- min(w[-1], na.rm = TRUE)
      ### Update the progress bar.
      if (w[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = w[1])
      }
      return(min.bcs)
    })
    # Max normalized BCS per signature.
    data[, 1] <- data[, 1] + n.rows
    max.bcscore <- apply(data, 1, function(x) {
      max.bcs <- max(x[-1], na.rm = TRUE)
      ### Update the progress bar.
      if (x[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = x[1])
      }
      return(max.bcs)
    })
    # NA proportion per signature.
    data[, 1] <- data[, 1] + n.rows
    prop.na <- apply(data, 1, function(y) {
      nas <- round(sum(is.na(y[-1]))/length(y[-1]), digits = 2)
      ### Update the progress bar.
      if (y[1]%%bins == 0) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = y[1])
      }
      return(nas)
    })
    # Residuals.
    normalized <- cbind(seq_len(n.rows) + (n.rows * i * 5) + (n.rows * (i - 1)),
                        bc@normalized[signatures, cells])
    # Create dataframe.
    stats <- data.frame(switch.point = switch.p, mean = mean.med.sd$mean,
                        median = mean.med.sd$median, sd = mean.med.sd$sd,
                        variance = variance.bcscore, min = min.bcscore,
                        max = max.bcscore, prop.na = prop.na,
                        row.names = signatures)
  } else {
    # Mean BCS per signature.
    mean.bc <- round(rowMeans(bc@data[signatures, cells], na.rm = TRUE),
                     digits = 2)
    # Residuals.
    normalized <- cbind(seq_len(n.rows) + (n.rows * (i - 1)),
                        bc@normalized[signatures, cells])
    # Create dataframe.
    stats <- data.frame(switch.point = switch.p, mean = mean.bc,
                        row.names = signatures)
  }
  # Residuals' mean.
  resid <- apply(normalized, 1, function(z) {
    res <- round(mean(z[-1], na.rm = TRUE), digits = 2)
    ### Update the progress bar.
    if (z[1]%%bins == 0 | z[1] == total) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, value = z[1])
    }
    return(res)
  })
  # Update dataframe.
  stats <- cbind(stats, data.frame(residuals.mean = resid, row.names = signatures))
  return(stats)
}

#' @title Returns the IDs that match the specified filter's values
#' @description This function subsets \code{df} to select only the entries that
#' match the specified \code{filter}'s \code{values} and returns the
#' corresponding IDs.
#' @name GetIDs
#' @param values User-supplied filtering vector.
#' @param filter Column name or number to subset by.
#' @param df \code{data.frame} with drug information. It must contain, at least,
#' two columns: \code{"IDs"} and \code{filter}.
#' @return A vector with the IDs that match the \code{filter}'s values.
#' @export

GetIDs <- function(values, filter, df = drugInfo) {
  # --- Checks ---
  # Check values.
  if (length(values) < 1 | !is.character(values)) {
    stop('values must be a character vector.')
  }
  # Check filter.
  if (length(filter) != 1) {
    stop('You must specify a single filter.')
  }
  if (is.character(filter) & !(filter %in% colnames(df))) {
    stop(paste('filter =', filter, 'is not a column of df.'))
  }
  if (is.numeric(filter) & (filter < 1 | filter > ncol(df))) {
    stop(paste('filter = ', filter, 'is out of range.'))
  }
  # Check df.
  if (class(df) != "data.frame") {
    stop('df must be a data.frame')
  }
  if (!("IDs" %in% colnames(df))) {
    stop('df must contain an "IDs" column.')
  }
  # --- Code ---
  upper.values <- toupper(values)
  selected <- subset(df, subset = toupper(df[[filter]]) %in% upper.values)
  if (filter == "drugs" & "preferred.drug.names" %in% colnames(df)) {
    synonyms <- subset(df, subset = toupper(df[["preferred.drug.names"]]) %in%
                         unique(toupper(selected[["preferred.drug.names"]])))
    selected <- unique(rbind(selected, synonyms))
  }
  ids <- unique(selected$IDs)
  not.found <- values[!(upper.values %in% toupper(df[[filter]]))]
  if (all(values %in% not.found)) {
    stop('No sig ID was found for any of the elements in values.')
  } else if (length(not.found) > 0) {
    filtername <- gsub(pattern = '"', replacement = '',
                       x = deparse(substitute(filter)))
    warning(paste0('sig IDs were not found for ', length(not.found), ' out of ',
                   length(values), " ", filtername, ': ',
                   paste0(not.found, collapse = ", "), "."))
  }
  return(ids)
}

#' @title Returns a dataframe with information about the input drugs
#' @description This function searches the input drugs in the pre-loaded
#' \code{beyondcell} matrices and returns a dataframe with drug information,
#' including drug synonyms and MoAs.
#' @name FindDrugs
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param x A character vector with drug names and/or sig IDs.
#' @param na.rm Logical. Should \code{x} entries with no available drug 
#' information be removed from the final output?
#' @details The output \code{data.frame} has the following columns:
#' \itemize{
#' \item{\code{original.names}}: Input drug names.
#' \item{\code{bc.names}}: Drug names used in \code{bc}.
#' \item{\code{preferred.drug.names}}: Standard drug names.
#' \item{\code{drugs}}: Other drug names.
#' \item{\code{IDs}}: Signature IDs.
#' \item{\code{preferred.and.sigs}}: \code{preferred.drug.names} (or alternatively
#' \code{bc.names}) and \code{IDs}. Used as title in \code{beyondcell} plots.
#' \item{\code{MoAs}}: Mechanism(s) of action.
#' }
#' @return A \code{data.frame}.
#' @examples
#' @export

FindDrugs <- function(bc, x, na.rm = TRUE) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check x.
  if (!is.character(x)) stop('x must be a character vector.')
  # Check na.rm.
  if (length(na.rm) != 1 | !is.logical(na.rm)) {
    stop('na.rm must be TRUE or FALSE.')
  }
  # --- Code ---
  # bc signatures.
  sigs <- rownames(bc@normalized)
  # Match x with bc signatures and get the indexes of matching elements.
  indexes <- lapply(x, function(y) {
    idx <- match(toupper(y), table = toupper(sigs), nomatch = 0)
    if (idx == 0) {
      idx <- unique(match(drugInfo$IDs[drugInfo$drugs == toupper(y)],
                          table = sigs))
    }
    return(idx[!is.na(idx)])
  })
  # Original names (x) and bc names (sigs).
  df <- data.frame(original.names = unlist(sapply(seq_along(x), function(i) {
    rep(x[i], times = length(indexes[[i]]))
  })), IDs = unlist(sapply(indexes, function(z) sigs[z])))
  df.not.found <- !(x %in% df$original.names)
  if (any(df.not.found)) {
    empty.df <- data.frame(original.names = x[df.not.found],
                           IDs = rep(NA, sum(df.not.found)))
    df <- rbind(df, empty.df)
  }
  # Get the names and pathways of the selected signatures.
  info <- subset(drugInfo, subset = IDs %in% df$IDs)
  if (all(dim(info) != 0)) {
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(w) {
      paste(na.omit(unique(w)), collapse = ", ")
    })
  }
  info.not.found <- !(df$IDs %in% drugInfo$IDs)
  if (any(info.not.found)) {
    empty.info <- matrix(rep(NA, times = sum(info.not.found)*6), ncol = 6,
                         dimnames = list(1:sum(info.not.found),
                                         colnames(drugInfo)))
    info <- rbind(info, as.data.frame(empty.info))
  }
  # Merge df and info.
  df <- unique(merge(df, info[, c("IDs", "drugs", "preferred.drug.names",
                                  "MoAs")], by = "IDs", all.x = TRUE))
  # Add bc.names column and remove names that are not sig IDs from sig_id
  # column.
  df$bc.names <- df$IDs
  df$IDs[!startsWith(df$IDs, prefix = "sig_")] <- NA
  # Create preferred.and.sigs column: Preferred name and sig_id.
  df$preferred.and.sigs <- sapply(1:nrow(df), function(j) {
    return(ifelse(test = !is.na(df$preferred.drug.names[j]),
                  yes = paste0(df$preferred.drug.names[j],
                               paste0(" (", df$IDs[j], ")")),
                  no = df$bc.names[j]))
  })
  # Reorder df.
  rows <- unlist(lapply(x, function(entry) which(df$original.names == entry)))
  cols <- c("original.names", "bc.names", "preferred.drug.names", "drugs", "IDs",
            "preferred.and.sigs", "MoAs")
  df <- df[rows, cols]
  # If na.rm = TRUE, remove rows with NAs in "preferred.drug.names" and "drugs" 
  # fields.
  if (na.rm) df <- df[rowSums(is.na(df[, 3:4])) < 2, ]
  return(df)
}
