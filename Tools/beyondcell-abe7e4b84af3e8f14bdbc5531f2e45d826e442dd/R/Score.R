#' @title Computes the BCS
#' @description This function computes the beyondcell score (BCS) and returns an
#' object of class \code{\link[beyondcell]{beyondcell}}.
#' @name bcScore
#' @import Seurat
#' @import scales
#' @param sc \code{\link[Seurat]{Seurat}} object or expression matrix.
#' @param gs \code{\link[beyondcell]{geneset}} object.
#' @param expr.thres Minimum fraction of signature genes that must be
#' expressed in a cell to compute its BCS. Cells with a number of expressed
#' genes below this fraction will have a \code{NaN} BCS.
#' @return A \code{beyondcell} object.
#' @examples
#' @export

bcScore <- function(sc, gs, expr.thres = 0.1) {
  # --- Checks ---
  # Check if sc is a Seurat object or an expression matrix.
  if ("Seurat"%in% class(sc)) {
    input <- "Seurat object"
    default <- Seurat::DefaultAssay(sc)
    if(default %in% c("RNA", "SCT")) {
      if("data" %in% slotNames(sc@assays[[default]])) {
        message(paste0('Using ', default, ' assay as input.'))
        expr.matrix <- as.matrix(Seurat::GetAssayData(sc, slot = "data",
                                                      assay = default))
      } else {
        stop('Default assay must include a normalized data (@data) slot.')
      }
    } else {
      stop('Seurat default assay must be either RNA or SCT.')
    }
  } else if ("matrix" %in% class(sc) & is.numeric(sc)) {
    input <- "expression matrix"
    warning(paste('Using count matrix as input. Please, check that this matrix',
                  'is normalized and unscaled.'))
    expr.matrix <- sc
    sc <- Seurat::CreateSeuratObject(expr.matrix)
  } else stop(paste('sc must be either a Seurat object or a single-cell',
                    'expression matrix.'))
  # Check if gs is a geneset object.
  if (class(gs) != "geneset") stop('gs must be a geneset object.')
  # Check expr.thres.
  if (length(expr.thres) != 1 | expr.thres[1] < 0 | expr.thres[1] > 1) {
    stop('expr.thres must be a positive number between 0 and 1.')
  }
  # Check that gene names are in the same format.
  sc.gene.case <- names(which.max(CaseFraction(rownames(expr.matrix))))
  gs.gene.case <- names(which.max(CaseFraction(unique(unlist(gs@genelist)))))
  if (sc.gene.case != gs.gene.case) {
    warning(paste0('gs genes are ', sc.gene.case, ' and sc genes are ',
                   gs.gene.case, '. Please check your ', input,
                   ' and translate the genes if necessary.'))
  }
  # --- Code ---
  # Create beyondcell object.
  bc <- beyondcell(expr.matrix = expr.matrix, meta.data = sc@meta.data,
                   SeuratInfo = list(reductions = sc@reductions),
                   regression = list(order = rep("", 2), vars = NULL,
                                     order.background = rep("", 2)),
                   n.genes = gs@n.genes, mode = gs@mode, thres = expr.thres)
  # Convert genes in expr.matrix and gs to lowercase.
  rownames(expr.matrix) <- tolower(rownames(expr.matrix))
  gs@genelist <- lapply(gs@genelist, function(x) {
    if (is.null(x[["up"]])) up <- NULL
    else up <- tolower(x[["up"]])
    if (is.null(x[["down"]])) down <- NULL
    else down <- tolower(x[["down"]])
    list(up = up, down = down)
  })
  # Genes in expr.matrix.
  genes <- rownames(expr.matrix)
  # Progress bar.
  len.gs <- length(gs@genelist)
  total <- len.gs + (length(gs@mode) * len.gs)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  bins <- ceiling(total / 100)
  # Below expr.thres (common for "up" and "down" modes).
  below.thres <- t(sapply(1:len.gs, function(i) {
    ### All genes in signature ("up" and "down").
    all.genes <- unique(unlist(gs@genelist[[i]]))
    ### Subset expr.matrix.
    sub.expr.matrix <- expr.matrix[genes %in% all.genes, , drop = FALSE]
    ### Number of genes expressed (expression > 0).
    n.expr.genes <- apply(sub.expr.matrix, 2, function(x) sum(x > 0))
    ### Update the progress bar.
    if (i%%bins == 0) {
      Sys.sleep(0.1)
      setTxtProgressBar(pb, value = i)
    }
    ### Is n.expr.genes < all.genes * expr.thres?
    return(n.expr.genes < (length(all.genes) * expr.thres))
  }))
  rownames(below.thres) <- names(gs@genelist)
  below.thres <- below.thres[, colnames(expr.matrix)]
  # If all cells are below the threshold, remove that signature and raise a
  # warning
  nan.rows.idx <- which(rowSums(below.thres) == ncol(below.thres))
  if (length(nan.rows.idx) > 0) {
    warning(paste0("The following signatures have no cells that pass the ", 
                   "expr.thres and will be removed: ", 
                   paste0(rownames(below.thres)[nan.rows.idx], collapse = ", "), 
                   "."))
    below.thres <- below.thres[-nan.rows.idx, ]
    gs@genelist <- gs@genelist[-nan.rows.idx]
  }
  # BCS.
  bcs <- lapply(seq_along(gs@mode), function(j) {
    score <- t(sapply(1:length(gs@genelist), function(k) {
      ### Common genes between the expr.marix and each signature.
      common.genes <- unique(intersect(genes, gs@genelist[[k]][[gs@mode[j]]]))
      ### Subset expr.matrix with common.genes.
      sub.expr.matrix <- expr.matrix[common.genes, , drop = FALSE]
      ### Sum expression of those genes for each cell.
      sum.expr <- colSums(sub.expr.matrix)
      ### Raw score (mean).
      raw <- colMeans(sub.expr.matrix)
      ### Stdev (for BCS normalization).
      sig.stdev <- apply(sub.expr.matrix, 2, sd)
      ### Normalized BCS.
      norm.score <- raw * ((sum.expr - sig.stdev)/(raw + sig.stdev))
      ### Update the progress bar.
      step <- len.gs + (j - 1) * len.gs + k + length(nan.rows.idx)
      if (step%%bins == 0 | step == total) {
        Sys.sleep(0.1)
        setTxtProgressBar(pb, value = step)
      }
      return(norm.score)
    }))
    rownames(score) <- names(gs@genelist)
    return(score[, colnames(expr.matrix), drop = FALSE])
  })
  names(bcs) <- gs@mode
  # Operate "up" and "down" scores.
  if (length(gs@mode) == 2) {
    ### Which cells are NaN in both modes.
    nan.cells <- is.na(bcs[["up"]]) & is.na(bcs[["down"]])
    ### Change NaNs to 0 so we can operate.
    bcs[["up"]][is.na(bcs[["up"]])] <- 0
    bcs[["down"]][is.na(bcs[["down"]])] <- 0
    ### Operate "up" and "down" scores.
    scoring.matrix <- bcs[["up"]] - bcs[["down"]]
    ### Restore NaNs.
    scoring.matrix[nan.cells | below.thres] <- NaN
  } else {
    scoring.matrix <- bcs[[gs@mode]]
    ### Restore NaNs.
    scoring.matrix[below.thres] <- NaN
    ### If mode is "down", values must be negative.
    if (gs@mode == "down") {
      scoring.matrix <- -1 * scoring.matrix[, , drop = FALSE]
    }
  }
  # If genesets were obtained in a TREATED vs CONTROL comparison, invert BCS
  # sign (exclude pathways).
  if (gs@comparison == "treated_vs_control") {
    not.paths <- which(!(rownames(scoring.matrix) %in% names(pathways)))
    paths <- which(rownames(scoring.matrix) %in% names(pathways))
    if (length(paths) > 0) {
      scoring.matrix <- rbind((scoring.matrix[not.paths, , drop = FALSE] * -1),
                              scoring.matrix[paths, , drop = FALSE])
    } else {
      scoring.matrix <- (-1) * scoring.matrix[, , drop = FALSE]
    }
  }
  slot(bc, "normalized") <- slot(bc, "data") <- round(scoring.matrix, digits = 2)
  # Scale the final matrix [0, 1] for each signature (for visualization).
  scaled.matrix <- t(apply(bc@normalized, 1, scales::rescale, to = c(0, 1)))
  slot(bc, "scaled") <- round(scaled.matrix, digits = 2)
  # Compute the switch point.
  switch.point <- SwitchPoint(bc)
  slot(bc, "switch.point") <- switch.point
  # Close the progress bar.
  close(pb)
  return(bc)
}

#' @title Returns the fraction of each case type in the input
#' @description This function computes the fraction of of each case type
#' (uppercase, lowercase or capitalized) in a character vector.
#' @name CaseFraction
#' @import useful
#' @param x Character vector.
#' @return A named numeric vector with the fraction of each case type.
#' @examples
#' @export

CaseFraction <- function(x) {
  # --- Checks ---
  if(!is.character(x)) stop('x must be a character vector.')
  # --- Code ---
  perc <- sapply(c("upper", "lower", "capitalized"), function(case) {
    if (case == "upper" | case == "lower") {
      p <- sum(useful::find.case(x, case = case))/length(x)
    } else {
      splitted <- strsplit(x, split = "")
      first.letter <- sapply(splitted, `[[`, 1)
      rest.letters <- sapply(splitted, function(y) paste0(y[2:length(y)],
                                                          collapse = ""))
      p <- sum(useful::upper.case(first.letter) &
                 useful::lower.case(rest.letters))/length(x)
    }
    return(round(p, digits = 2))
  })
  names(perc)[1:2] <- paste0("in ", names(perc)[1:2], "case")
  return(perc)
}

#' @title Computes the switch point
#' @description This function computes the switch point of the signatures of a
#' given \code{\link[beyondcell]{beyondcell}} object. The switch point is the
#' (subsetted and/or regressed) scaled beyondcell score (BCS) that corresponds
#' to the point in which the normalized BCS in \code{beyondcell@@data} switchs
#' from negative (insensitive) to positive (sensitive) values. The closer to 0,
#' the more sensitive are the cells to a given drug.
#' @name SwitchPoint
#' @param bc \code{beyondcell} object.
#' @return A named vector with the swich points of the signatures in \code{bc}.
#' @examples
#' @export

SwitchPoint <- function(bc) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check bc@mode.
  if (any(!(bc@mode %in% c("up", "down"))) | all(is.null(bc@mode))) {
    stop('Incorrect mode.')
  }
  bc@mode <- unique(bc@mode)
  # --- Code ---
  # Signatures in bc.
  sigs <- rownames(bc@normalized)
  # Cells in bc.
  cells <- colnames(bc@normalized)
  # If bc@mode == "up", all normalized BCS will be positive and all switch
  # points will be 0.
  if (all(bc@mode == "up")) {
    switch.point <- rep(0, times = length(sigs))
    # If bc@mode == "down", all normalized BCS will be negative and all switch
    # points will be 1.
  } else if (all(bc@mode == "down")) {
    switch.point <- rep(1, times = length(sigs))
    # If bc@mode == c("up", "down")...
  } else if(length(bc@mode) == 2) {
    ### Create a list with an entry for each signature. If the entry has
    ### length == 1, it is the final switch point. If it has length == 2, it
    ### corresponds to the indexes of the normalized BCS closest to 0.
    indexes <- lapply(sigs, FUN = function(x) {
      ### Subset the normalized BCS in @data.
      m <- bc@data[x, cells]
      ### If all values are NaN, return NaN.
      if (all(is.na(m))) return(NaN)
      ### Else, remove NaN values.
      m.nona <- na.omit(m)
      ### If all normalized BCS are positive, the switch point is 0.
      if (all(m.nona >= 0)) return(0)
      ### If all normalized BCS are negative, the switch point is 1.
      else if (all(m.nona <= 0)) return(1)
      ### If there are positive and negative normalized BCS...
      else {
        exact.0 <- m.nona == 0
        ### If any of the normalized BCS == 0, return its index (we repeat
        ### it twice to indicate it's an index and not the final switch point).
        if (any(exact.0)) return(rep(which(m == 0)[1], times = 2))
        ### Else, get the indexes of the values closer to 0.
        else {
          lower.bound <- which(m == max(m.nona[m.nona <= 0]))[1]
          upper.bound <- which(m == min(m.nona[m.nona >= 0]))[1]
          return(c(lower.bound, upper.bound))
        }
      }
    })
    names(indexes) <- sigs
    ### For each indexes entry, return a vector with the switch point.
    switch.point <- sapply(names(indexes), function(y) {
      ### If length(entry) == 2, the values are indexes. We subset those
      ### indexes in the scaled bcsores and take the arithmetic mean.
      if (length(indexes[[y]]) == 2) {
        return(round(sum(bc@scaled[y, indexes[[y]]])/2, digits = 2))
        ### If length(entry) == 1, the value is the final switch point.
      } else {
        return(indexes[[y]])
      }
    })
  }
  names(switch.point) <- sigs
  return(switch.point)
}
