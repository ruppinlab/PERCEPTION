#' @title Plots a UMAP reduction coloured by metadata information
#' @description This function returns a {\link[ggplot2]{ggplot2}} object with
#' a UMAP reduction coloured by the specified metadata column.
#' @name bcClusters
#' @import Seurat
#' @import ggplot2
#' @import scales
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param idents Name of the metadata column to colour by.
#' @param UMAP UMAP reduction to plot. Either \code{"beyondcell"}, computed
#' using \code{\link[beyondcell]{bcUMAP}}, or \code{"Seurat"}, obtained using
#' \code{Seurat}'s functions.
#' @param factor.col Logical indicating if \code{idents} column is a factor or
#' not. Set \code{factor.col = FALSE} if \code{idents} is a numeric column (such
#' as \code{percent.mt} or \code{nFeature_RNA}).
#' @param ... Other arguments passed to \code{\link[Seurat]{DimPlot}}.
#' @return A \code{ggplot2} object with the UMAP reduction coloured by
#' \code{idents}.
#' @examples
#' @export

bcClusters <- function(bc, idents, UMAP = "beyondcell", factor.col = TRUE,
                       ...) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (length(idents) != 1) {
    stop('Idents must be a single metadata column.')
  }
  if (idents %in% colnames(bc@meta.data)) {
    meta <- bc@meta.data[colnames(bc@scaled), idents, drop = FALSE]
  } else {
    stop('Idents not found.')
  }
  # Check UMAP.
  if (UMAP == "beyondcell") {
    if (length(bc@reductions) == 0) {
      stop('You must precompute beyondcell\'s UMAP projection using bcUMAP().')
    }
    sc <- Seurat::CreateSeuratObject(bc@scaled)
    reduction <- bc@reductions
  } else if (UMAP == "Seurat") {
    if (length(bc@SeuratInfo$reductions) == 0) {
      stop('No UMAP projection available for your Seurat\'s object.')
    }
    sc <- Seurat::CreateSeuratObject(bc@expr.matrix[, colnames(bc@scaled)])
    reduction <- bc@SeuratInfo$reductions
  } else {
    stop('Incorrect UMAP argument. Please use either "Seurat" or "beyondcell".')
  }
  # Check factor.col.
  if (length(factor.col) != 1 | !is.logical(factor.col)) {
    stop('factor.col must be TRUE or FALSE.')
  }
  # --- Code ---
  # Add metadata.
  sc <- Seurat::AddMetaData(sc, metadata = meta)
  # Add reductions.
  sc@reductions <- reduction
  # Plot.
  if (factor.col) {
    # Set Idents.
    Seurat::Idents(sc) <- idents
    p <- Seurat::DimPlot(sc, reduction = "umap", ...) + ggplot2::theme_minimal()
  } else {
    p <- Seurat::FeaturePlot(sc, reduction = "umap", features = idents, ...) +
      ggplot2::theme_minimal() + ggplot2::labs(title = NULL)
  }
  return(p)
}

#' @title Plots a histogram with the BCS of the signature of interest
#' @description This function drawns a histogram of beyondcell scores (BCS) for
#' each signature of interest. The plot can be a single histogram or a histogram
#' for each level found in \code{idents}.
#' @name bcHistogram
#' @import ggplot2
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param signatures Vector with the names of the signatures of interest. If
#' \code{signatures = "all"}, all signatures are selected.
#' @param idents Name of the metadata column of interest. If
#' \code{idents = NULL}, a single histogram with all BCS is drawn. On the other
#' hand, if \code{idents != NULL} a histogram for each level found in
#' \code{idents} will be drawn.
#' @return A list of \code{\link[ggplot2]{ggplot2}} histograms, one for each
#' signature of interest. In each histogram, the median, mean and sd are
#' reported. Also, the mean is indicated with a black dashed line and the median
#' with a red dashed line.
#' @examples
#' @export

bcHistogram <- function(bc, signatures, idents = NULL) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check signatures.
  if (!is.character(signatures)) {
    stop("Signatures must be a character vector.")
  }
  if (length(signatures) == 1 & signatures[1] == "all") {
    signatures <- rownames(bc@normalized)
    in.signatures <- rep(TRUE, times = nrow(bc@normalized))
  } else {
    in.signatures <- !is.null(signatures) &
      signatures %in% rownames(bc@normalized)
    if (all(!in.signatures)) {
      stop('None of the specified signatures were found.')
    } else if (any(!in.signatures)) {
      warning(paste0('These signatures were not found in bc: ',
                     paste0(signatures[!in.signatures], collapse = ", "), '.'))
    }
  }
  # Check idents.
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop("Idents must be a single metadata column.")
    }
    if (idents %in% colnames(bc@meta.data)) {
      meta <- paste(idents, "=", bc@meta.data[colnames(bc@normalized), idents,
                                              drop = TRUE])
    } else {
      stop('Idents not found.')
    }
  } else {
    ### If idents = NULL, all cells have the same metadata.
    meta <- rep("", times = ncol(bc@normalized))
  }
  # --- Code ---
  # Metadata levels.
  lvls <- levels(as.factor(meta))
  # Subset bc to the selected signatures.
  sub.bc <- bc@data[signatures[in.signatures], , drop = FALSE]
  # Get maximum and minimum normalized BCS (for common x axis in all plots).
  limits <- c(min(as.vector(sub.bc), na.rm = TRUE),
              max(as.vector(sub.bc), na.rm = TRUE))
  # Get info about drugs (their corresponding name in bc, the preferred name
  # used by beyondcell and the MoA).
  info <- FindDrugs(bc, x = signatures[in.signatures])
  # For each signature.
  p <- lapply(signatures[in.signatures], function(x) {
    ### Data frame of normalized BCS and metadata.
    sub.df <- na.omit(data.frame(bcscore = sub.bc[x, ],
                                 condition = as.factor(meta),
                                 row.names = colnames(sub.bc)))
    ### Statistics by metadata (mean, median and sd).
    stats.cond <- sapply(lvls, function(y) {
      stats <- round(Mean.Med.SD(subset(sub.df$bcscore, subset = sub.df$condition == y)),
                     digits = 2)
      return(stats)
    })
    stats.labels <- data.frame(label = apply(stats.cond, 2, function(z) {
      paste0(rownames(stats.cond), " = ", z, collapse = "\n")
    }), mean = stats.cond["mean", ], median = stats.cond["median", ],
    condition = colnames(stats.cond))
    ### Drug name and MoA
    if (x %in% info$IDs) {
      drug.and.MoA <- info[which(info$IDs == x), c("preferred.and.sigs", "MoAs")]
      drug.and.MoA[2] <- ifelse(test = drug.and.MoA[2] == "NA", yes = "",
                                no = drug.and.MoA[2])
    } else {
      drug.and.MoA <- c(x, "")
    }
    ### Plot
    hist <- ggplot(sub.df, aes(x = bcscore, fill = condition)) +
      geom_histogram(data = transform(sub.df, condition = NULL),
                     fill = "grey85", alpha = 0.3, binwidth = 1) +
      facet_wrap(vars(condition), ncol = 1) +
      geom_histogram(color = "#F0F0F0", alpha = 0.6, position = "identity",
                     binwidth = 1) + theme_minimal() +
      theme(legend.position = "none") +
      geom_label(data = stats.labels, hjust = 1, vjust = 1, size = 3,
                 x = Inf, y = Inf, label = stats.labels$label,
                 fill = scales::hue_pal()(length(lvls))) +
      labs(title = drug.and.MoA[1], subtitle = drug.and.MoA[2]) +
      geom_vline(data = stats.labels, aes(xintercept = mean),
                 linetype = "dashed") +
      geom_vline(data = stats.labels, aes(xintercept = median, color = "red"),
                 linetype = "dashed")
    return(hist)
  })
  return(p)
}

#' @title Plots a UMAP reduction coloured by BCS or gene expression values
#' @description This function returns a list of
#' \code{\link[patchwork]{patchwork}} or \code{\link[ggplot2]{ggplot2}} objects
#' with the desired UMAP reduction coloured by beyondcell scores (BCS) or gene
#' expression values.
#' @name bcSignatures
#' @import Seurat
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param UMAP UMAP reduction to plot. Either \code{"beyondcell"}, computed
#' using \code{\link[beyondcell]{bcUMAP}}, or \code{"Seurat"}, obtained using
#' \code{Seurat}'s functions.
#' @param signatures List with plot parameters to colour the UMAP by BCS:
#' \itemize{
#' \item{\code{values}:} {Vector with the names of the signatures of interest.
#' If \code{signatures[["values"]] = "all"}, all signatures are selected.}
#' \item{\code{colorscale}:} {Either a \code{viridis}, \code{RColorBrewer} or a
#' custom palette of 3 colours (low, medium and high) to colour all signatures'
#' plots. If \code{colorscale = NULL} (default), the plots are coloured using
#' \code{beyondcell}'s own palette.}
#' \item{\code{alpha}:} {Transparency level between 1 (not transparent) and 0
#' (fully transparent).}
#' \item{\code{na.value}:} {Colour to use for missing values (\code{NA}s).}
#' \item{\code{limits}:} {Vector with the desired limits for all signatures'
#' plots.}
#' \item{\code{center}:} {A single number indicating the centre of the
#' \code{colorscale} for all signatures' plots. If \code{center = NULL}
#' (default), the \code{center} for each signature is its switch point.}
#' \item{\code{breaks}:} {A single number indicating the break size of the
#' \code{colorscale}. Alternatively, it can be a vector with the desired breaks
#' (which don't have to be symmetric or equally distributed).}}
#' @param genes List with plot parameters to colour the UMAP by gene expression
#' values:
#' \itemize{
#' \item{\code{values}:} {Vector with the names of the genes of interest. If
#' \code{genes[["values"]] = "all"}, all genes are selected.}
#' \item{\code{limits}:} {Vector with the desired limits for all genes' plots.
#' If \code{limits = c(NA, NA)} (default), the \code{limits} are computed
#' automatically. See Details for  more information.}
#' \item{\code{share.limits}:} {Logical argument. If \code{share.limits = TRUE},
#' all genes' plots will have the same \code{limits}. If
#' \code{share.limits = FALSE} (default), each gene plot will have its own
#' \code{limits}. See Details for more information.}}
#' @param merged If \code{merged != NULL}, two signatures will be superposed in
#' the same plot. If \code{merged = "direct"}, the signatures are assumed to
#' have a direct relationship and the BCS will be added. On the other hand, if
#' \code{merged = "indirect"}, the signatures are assumed to have an indirect
#' relationship and their BCS will be substracted.
#' @param blend (From \code{\link[Seurat]{FeaturePlot}}) Scale and blend
#' expression values to visualise co-expression of two genes.
#' @param mfrow Numeric vector of the form \code{c(nr, nc)}. \code{nr}
#' corresponds to the number of rows and \code{nc} to the number of columns of
#' the grid in which the plots will be drawn. If you want to draw the plots
#' individually, set \code{mfrow = c(1, 1)}.
#' @param ... Other arguments passed to \code{FeaturePlot}.
#' @details When \code{genes[["limits"]] = c(NA, NA)}, \code{bcSignatures}
#' computes the limits automatically. You can make all plots share the same
#' limits by specifying \code{genes[["share.limits"]] = TRUE}, or make the
#' function to compute the limits individually for each gene with
#' \code{genes[["share.limits"]] = FALSE}. Moreover, if you specify a value for
#' \code{genes[["limits"]]}, \code{genes[["share.limits"]]} will automatically
#' set to \code{TRUE} and all plots will share those limits.
#' @return A list of \code{patchwork} (if \code{mfrow != c(1, 1)}) or
#' \code{ggplot2} objects (if \code{mfrow = c(1, 1)}) of the desired UMAP
#' reduction coloured by the BCS (for signatures) or gene expression values (for
#' genes).
#' @examples
#' @export

bcSignatures <- function(bc, UMAP = "beyondcell",
                         signatures = list(values = NULL, colorscale = NULL,
                                           alpha = 0.7, na.value = "grey50",
                                           limits = c(0, 1), center = NULL,
                                           breaks = 0.1),
                         genes = list(values = NULL, limits = c(NA, NA),
                                      share.limits = FALSE),
                         merged = NULL, blend = FALSE, mfrow = c(1, 1), ...) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check UMAP.
  if (UMAP == "beyondcell") {
    if (length(bc@reductions) == 0) {
      stop('You must precompute beyondcell\'s UMAP projection using bcUMAP().')
    }
    reduction <- bc@reductions
    cells <- subset(rownames(bc@meta.data),
                    subset = rownames(bc@meta.data) %in% colnames(bc@normalized))
  } else if (UMAP == "Seurat") {
    if (length(bc@SeuratInfo$reductions) == 0) {
      stop('No UMAP projection available for your Seurat\'s object.')
    }
    reduction <- bc@SeuratInfo$reductions
    cells <- rownames(bc@meta.data)
  } else {
    stop('Incorrect UMAP argument. Please use either "Seurat" or "beyondcell".')
  }
  # Check signatures' list values.
  default.sigs <- list(values = NULL, colorscale = NULL, alpha = 0.7,
                       na.value = "grey", limits = c(0, 1), center = NULL,
                       breaks = 0.1)
  selected.sigs <- names(signatures) %in% names(default.sigs)
  if (any(!selected.sigs)) {
    warning(paste0('Incorrect entries in signatures: ',
                   paste0(names(signatures)[!selected.sigs], collapse = ", ")))
  }
  signatures <- c(signatures,
                  default.sigs[!(names(default.sigs) %in% names(signatures))])
  # Check genes' list values.
  default.genes <- list(values = NULL, limits = c(NA, NA), share.limits = FALSE)
  selected.genes <- names(genes) %in% names(default.genes)
  if (any(!selected.genes)) {
    warning(paste0('Incorrect entries in genes: ',
                   paste0(names(genes)[!selected.genes], collapse = ", ")))
  }
  genes <- c(genes, default.genes[!(names(default.genes) %in% names(genes))])
  # Check signatures and genes' values (features).
  if (is.null(signatures[["values"]]) & is.null(genes[["values"]])) {
    stop('You must specify the signatures and/or genes of interest.')
  }
  # Check signatures' values.
  if (!is.null(signatures[["values"]])) {
    if (length(signatures[["values"]]) == 1 & signatures[["values"]][1] == "all") {
      signatures[["values"]] <- rownames(bc@normalized)
      in.signatures <- rep(TRUE, times = nrow(bc@normalized))
    } else {
      in.signatures <- signatures[["values"]] %in% rownames(bc@normalized)
      if (all(!in.signatures)) {
        stop('None of the specified signatures were found.')
      } else if (any(!in.signatures)) {
        warning(paste0('These signatures were not found in bc: ',
                       paste0(signatures[["values"]][!in.signatures],
                              collapse = ", "), '.'))
      }
    }
  } else {
    in.signatures <- NULL
  }
  # Check genes' values.
  if (!is.null(genes[["values"]])) {
    if (length(genes[["values"]]) == 1 & genes[["values"]][1] == "all") {
      genes[["values"]] <- rownames(bc@expr.matrix)
      in.genes <- rep(TRUE, times = nrow(bc@expr.matrix))
    } else {
      in.genes <- toupper(genes[["values"]]) %in% toupper(rownames(bc@expr.matrix))
      if (all(!in.genes)) {
        stop('None of the specified genes were found.')
      } else if (any(!in.genes)) {
        warning(paste0('These genes were not found in bc@expr.matrix: ',
                       paste0(genes[["values"]][!in.genes], collapse = ", "), '.'))
      }
    }
  } else {
    in.genes <- NULL
  }
  # Sigs, gene and features.
  sigs <- unique(signatures[["values"]][in.signatures])
  gene <- unique(genes[["values"]][in.genes])
  features <- c(sigs, gene)
  # Check signature's colorscale.
  signatures[["colorscale"]] <- get_colour_steps(signatures[["colorscale"]])
  # Check signatures' alpha, na.value and breaks -> inside center_scale_colour_stepsn().
  # Check signatures' limits.
  if (length(signatures[["limits"]]) != 2) {
    stop('Signatures\' limits must be a vector of length 2.')
  }
  if (!is.numeric(signatures[["limits"]]) | any(signatures[["limits"]] < 0)) {
    stop('Signatures\' limits must be numeric (>= 0).')
  }
  if (signatures[["limits"]][2] < signatures[["limits"]][1]) {
    warning(paste('Signatures\' upper limit is smaller than lower limit.',
                  'Sorting limits in increasing order.'))
    signatures[["limits"]] <- sort(signatures[["limits"]], decreasing = FALSE)
  }
  # Check signature's center.
  if (!is.null(signatures[["center"]])) {
    if (length(signatures[["center"]]) != 1 |
        !is.numeric(signatures[["center"]])) {
      stop('Signatures\' center must be a single number or NULL.')
    }
  }
  # Check genes' limits.
  if (length(genes[["limits"]]) != 2) {
    stop('Genes\' limits must be a vector of length 2.')
  }
  na.limits.genes <- is.na(genes[["limits"]])
  if (length(genes[["limits"]][!na.limits.genes]) > 0 &
      (!is.numeric(genes[["limits"]][!na.limits.genes]) |
       any(genes[["limits"]][!na.limits.genes] < 0))) {
    stop('Genes\' limits must be numeric (>= 0) or NAs.')
  }
  if (all(!na.limits.genes) & genes[["limits"]][2] < genes[["limits"]][1]) {
    warning(paste('Genes\' upper limit is smaller than lower limit.',
                  'Sorting limits in increasing order.'))
    genes[["limits"]] <- sort(genes[["limits"]], decreasing = FALSE)
  }
  # Check genes' share.limits.
  if (length(genes[["share.limits"]]) != 1 | !is.logical(genes[["share.limits"]])) {
    stop('genes$share.limits must be TRUE or FALSE.')
  }
  if (!genes[["share.limits"]] &
      !identical(genes[["limits"]], c(NA, NA))) {
    warning(paste('Genes\' limits were specified, setting',
                  'genes[["share.limits"]] = TRUE.'))
  }
  # Check merged.
  if (!is.null(merged)) {
    if (length(merged) != 1 | !(merged %in% c("direct", "indirect"))) {
      stop(paste('Incorrect merged value: It must be either NULL, "direct"',
                 'or "indirect".'))
    }
    if (length(features) != 2 | length(signatures[["values"]]) != 2) {
      stop('When merged != NULL, the number of signatures must be exactly 2.')
    } else if (any(!(features %in% sigs))) {
      stop(paste('The merged features must be signatures. For blending genes,',
                 'please use blend = TRUE.'))
    }
    merged.symbol <- ifelse(test = merged == "direct", yes = " + ", no = " - ")
    merged.sigs <- paste0(sigs, collapse = merged.symbol)
    rm.NA <- FALSE
  } else {
    merged.symbol <- ""
    merged.sigs <- sigs
    rm.NA <- TRUE
  }
  # Check blend.
  if (length(blend) != 1 | !is.logical(blend)) {
    stop('blend must be TRUE or FALSE.')
  }
  if (blend) {
    if (length(features) != 2 | length(genes[["values"]]) != 2) {
      stop('When blend = TRUE, the number of genes must be exactly 2.')
    } else if (any(!(features %in% gene))) {
      stop(paste('The blended features must be genes. For merging signatures,',
                 'please use merged argument.'))
    }
  }
  # Check mfrow.
  if (length(mfrow) != 2 | any(mfrow < 0) | any(mfrow%%1 != 0)) {
    stop('mfrow must be a vector of two integers > 0.')
  }
  # --- Code ---
  # If blend = TRUE, plot a Seurat::FeaturePlot.
  if (blend) {
    ### Create fake seurat object.
    sc <- Seurat::CreateSeuratObject(bc@expr.matrix)
    ### Add reductions.
    sc@reductions <- reduction
    ### Plot.
    p <- suppressMessages(
      Seurat::FeaturePlot(sc, features = gene, blend = TRUE, combine = FALSE))
    # Else...
  } else {
    # Get info about drugs (their corresponding name in bc, the preferred name
    # used by beyondcell and the MoA).
    info <- tryCatch(suppressWarnings(FindDrugs(bc, x = sigs, na.rm = rm.NA)),
                     error = function(cond) data.frame())
    ### If we want to merge signatures, we must recompute the bc object using
    ### the added or substracted bc@normalized BCS.
    if (!is.null(merged)) {
      if (merged == "direct") {
        merged.bcscores <- colSums(bc@normalized[sigs, cells, drop = FALSE],
                                   na.rm = TRUE)
      } else if (merged == "indirect") {
        merged.bcscores <- colMinus(bc@normalized[sigs, cells, drop = FALSE],
                                    na.rm = TRUE)
      }
      bc@data <- matrix(merged.bcscores, nrow = 1, dimnames = list(merged.sigs, cells))
      bc <- suppressMessages(bcRecompute(bc, slot = "data"))
      features <- merged.sigs
    }
    ### Join scaled BCS and gene expression values for selected features.
    full.matrix <- rbind(bc@scaled[merged.sigs, cells, drop = FALSE],
                         bc@expr.matrix[gene, cells,
                                        drop = FALSE])[features, , drop = FALSE]
    ### Signature's center. If center == NULL, set center to switch.points.
    if (is.null(signatures[["center"]])) {
      center.sigs <- bc@switch.point[merged.sigs]
    } else {
      center.sigs <- setNames(rep(signatures[["center"]],
                                  times = length(merged.sigs)), merged.sigs)
    }
    ### Gene's colours.
    if ((genes[["share.limits"]] | !identical(genes[["limits"]], c(NA, NA))) &
        !is.null(genes[["values"]])) {
      if (any(na.limits.genes)) {
        if (!is.null(merged)) range.genes <- pretty(full.matrix)
        else range.genes <- pretty(full.matrix[gene, ])
        genes[["limits"]][na.limits.genes] <- c(min(range.genes),
                                                max(range.genes))[na.limits.genes]
      }
      colors.genes <- ggplot2::scale_colour_gradient(low = "lightgrey", high = "blue",
                                                     limits = genes[["limits"]])
    } else {
      colors.genes <- NULL
    }
    ### Create fake seurat object.
    sc <- Seurat::CreateSeuratObject(full.matrix)
    ### Add reductions.
    sc@reductions <- reduction
    ### Plot for each signtature/gene...
    p <- lapply(features, function(y) {
      ### Merged.
      if (!is.null(merged)) {
        ids <- unlist(strsplit(y, split = merged.symbol, fixed = TRUE))
      } else ids <- y
      ### Drug name and MoA.
      if (any(ids %in% info$IDs)) {
        drug.and.MoA <- info[which(info$IDs %in% ids),
                             c("preferred.and.sigs", "MoAs")]
        ### When merged != NULL...
        if (nrow(drug.and.MoA) > 1) {
          ### Paste both drug names and MoAs. If MoAs are the same, just print
          ### them one time
          drug.and.MoA <- as.data.frame(t(apply(drug.and.MoA, 2, FUN = function(z) {
            paste0(unique(z), collapse = merged.symbol)
          })))
        }
        drug.and.MoA[, 2] <- BreakString(drug.and.MoA[, 2]) ### Format subtitle.
      } else {
        ### Empty subtitle.
        drug.and.MoA <- c(paste0(ids, collapse = merged.symbol), "")
      }
      ### Colours (depending on wether y is a signature or a gene).
      if (all(ids %in% sigs)) {
        ### Binned colorscale centred around the switch point or the specified
        ### center.
        colors <- center_scale_colour_stepsn(full.matrix[y, ],
                                             colorscale = signatures[["colorscale"]],
                                             alpha = signatures[["alpha"]],
                                             na.value = signatures[["na.value"]],
                                             limits = signatures[["limits"]],
                                             center = center.sigs[y],
                                             breaks = signatures[["breaks"]])
      } else {
        ### Continuous colorscale with default Seurat colours.
        colors <- colors.genes
      }
      ### Plot.
      fp <- suppressMessages(
        Seurat::FeaturePlot(sc, features = gsub(pattern = "_", replacement = "-",
                                                x = y), combine = FALSE, ...)[[1]] +
          colors + ggplot2::labs(title = drug.and.MoA[1], subtitle = drug.and.MoA[2]))
      return(fp)
    })
  }
  # If mfrow = c(1, 1), return a list ggplots.
  if (identical(mfrow, c(1, 1))) {
    final.p <- p
    # Else, wrap plots according to mfrow and return a list of patchworks.
  } else {
    ncol.nrow <- mfrow[1] * mfrow[2]
    final.p <- lapply(1:ceiling(length(p)/ncol.nrow), function(i) {
      ### Range.
      start <- ((i - 1) * ncol.nrow) + 1
      end <- min(i * ncol.nrow, length(p))
      sub.p <- patchwork::wrap_plots(p[start:end], nrow = mfrow[1],
                                     ncol = mfrow[2])
      return(sub.p)
    })
  }
  return(final.p)
}

#' @title Plots a violindot of the BCS grouped by cell cycle phase
#' @description This function drawns, for each signature of interest, a
#' \code{\link[see]{geom_violindot}} plot of the beyondcell scores (BCS) grouped
#' by cell cycle phase (G1, G2M or S). Note that this information must be
#' present in \code{bc@@meta.data} and can be obtained using
#' \code{\link[Seurat]{CellCycleScoring}}.
#' @name bcCellCycle
#' @import ggplot2
#' @importFrom see geom_violindot
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param signatures Vector with the names of the signatures of interest. If
#' \code{signatures = "all"}, all signatures are selected.
#' @return A list of \code{\link[ggplot2]{ggplot2}} violindots, one for each
#' signature of interest. In each violindot, the BCS are grouped in G1, G2M or S
#' phase groups.
#' @examples
#' @export

bcCellCycle <- function(bc, signatures) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check that "Phase" column is present in bc@meta.data.
  if (!("Phase" %in% colnames(bc@meta.data))) {
    stop("Cell cycle information not present in bc@meta.data.")
  }
  # Check signatures.
  if (!is.character(signatures)) {
    stop("Signatures must be a character vector.")
  }
  if (length(signatures) == 1 & signatures[1] == "all") {
    signatures <- rownames(bc@normalized)
    in.signatures <- rep(TRUE, times = nrow(bc@normalized))
  } else {
    in.signatures <- !is.null(signatures) & signatures %in% rownames(bc@normalized)
    if (all(!in.signatures)) {
      stop('None of the specified signatures were found.')
    } else if (any(!in.signatures)) {
      warning(paste0('These signatures were not found in bc: ',
                     paste0(signatures[!in.signatures], collapse = ", "), '.'))
    }
  }
  # --- Code ---
  # Cells in beyondcell object.
  cells <- subset(rownames(bc@meta.data),
                  subset = rownames(bc@meta.data) %in% colnames(bc@normalized))
  # Get info about drugs (their corresponding name in bc, the preferred name
  # used by beyondcell and the MoA).
  info <- FindDrugs(bc, x = signatures[in.signatures])
  # For each signature...
  p <- lapply(signatures[in.signatures], function(x) {
    ### Data frame of normalized BCS and phase metadata.
    sub.df <- na.omit(data.frame(bcscore = bc@data[x, cells],
                                 phase = bc@meta.data[cells, "Phase"],
                                 row.names = cells))
    ### Drug name and MoA.
    if (x %in% info$IDs) {
      drug.and.MoA <- info[which(info$IDs == x),
                           c("preferred.and.sigs", "MoAs")]
      drug.and.MoA[, 2] <- BreakString(drug.and.MoA[, 2]) ### Format subtitle.
    } else {
      drug.and.MoA <- c(x, "")
    }
    ### Plot.
    violindot <- ggplot(sub.df, aes(x = phase, y = bcscore, fill = phase)) +
      see::geom_violindot() + theme_minimal() + theme(legend.position = "none") +
      labs(title = drug.and.MoA[1], subtitle = drug.and.MoA[2])
    return(violindot)
  })
  return(p)
}

#' @title Drawns a 4 squares plot
#' @description This function drawns a 4 square plot of the drug signatures
#' present in a \code{\link[beyondcell]{beyondcell}} object. A 4 squares plot
#' consists in a scatter plot of the residuals' means (x axis) vs the switch
#' points (y axis) of a specified cluster (either a therapeutic cluster or a
#' group defined by experimental condition or phenotype). 4 quadrants are
#' highlighted: the top-left and bottom-right corners contain the drugs to which
#' all selected cells are least/most sensistive, respectively. The centre
#' quadrants show the drugs to which these cells are differentially insensitive
#' or sensitive when compared to the other clusters.
#'
#' x cut-offs: first and last deciles; y cut-offs: 0.1, 0.4, 0.6 and 0.9.
#' @name bc4Squares
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#' @param bc \code{beyondcell} object.
#' @param idents Name of the metadata column of interest.
#' @param lvl Character vector with the \code{idents}' level(s) of interest. If
#' \code{lvl = NULL}, all levels will be plotted.
#' @param top Number of top drugs per quadrant to be labelled.
#' @param topnames Character vector with additional interesting drugs or 
#' pathways to be labeled (either their names or sig IDs).
#' @param force (From \code{\link[ggrepel]{ggrepel}}) Force of repulsion between
#' overlapping text labels. Defaults to 1.
#' @param alpha Transparency level between 1 (not transparent) and 0 (fully
#' transparent).
#' @param pt.size Point size.
#' @param ... Other arguments passed to \code{\link[ggrepel]{geom_text_repel}}.
#' @details This function returns a list of \code{\link[ggplot2]{ggplot2}}
#' objects, one per each \code{lvl}. Note that residuals' means are different
#' for each level while swicth points are signature-specific. So, x axis will
#' vary and y axis will remain constant accross all plots.
#' @return A list of \code{ggplot2} objects, one per \code{lvl}.
#' @examples
#' @export

bc4Squares <- function(bc, idents, lvl = NULL, top = 3, topnames = NULL, 
                       force = 1, alpha = 0.7, pt.size = 3, ...) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check idents.
  if (length(idents) != 1) stop('Idents must be a single metadata column.')
  if (!(idents %in% names(bc@ranks))) {
    stop(paste0('$', idents, ' not found in bc@ranks.'))
  } else if (idents == "general") {
    stop('General rank can\'t be used in bc4Squares(). All residuals are 0.')
  }
  # Check lvl.
  if (is.null(lvl)) {
    lvl <- unique(bc@meta.data[, idents])
    in.lvl <- rep(TRUE, times = length(lvl))
  } else {
    in.lvl <- lvl %in% unique(bc@meta.data[, idents])
    if (all(!in.lvl)) {
      stop(paste0('None of the specified levels were found in ', idents, '.'))
    } else if (any(!in.lvl)) {
      warning(paste0('The following levels were not found in ', idents, ': ',
                     paste0(lvl[!in.lvl], collapse = ", "), '.'))
    }
  }
  # Check top.
  if (length(top) != 1 | top[1]%%1 != 0 | top[1] < 0) {
    stop('top must be a single integer >= 0.')
  }
  # Check topnames.
  if (!is.null(topnames)) {
    in.topnames <- toupper(topnames) %in% drugInfo$drugs |
      tolower(topnames) %in% drugInfo$IDs |
      toupper(topnames) %in% toupper(rownames(bc@normalized))
    if (all(!in.topnames)) {
      warning('None of the specified topname drugs were found in bc.')
    } else if (any(!in.topnames)) {
      warning(paste0('The following topname drugs were not found in bc: ' ,
                     paste0(topnames[!in.topnames], collapse = ", "), '.'))
    }
  } else {
    in.topnames <- NULL
  }
  # Check force.
  if (length(force) != 1 | force[1] < 0) {
    stop('force must be a single number >= 0.')
  }
  # Check alpha.
  if (length(alpha) != 1 | alpha[1] < 0 | alpha[1] > 1) {
    stop('alpha must be a single number between 0 and 1.')
  }
  # Check pt.size.
  if (length(pt.size) != 1 | !is.numeric(pt.size)) {
    stop('pt.size must be a single number.')
  }
  # --- Code ---
  # Get info about drugs (their corresponding name in bc, the preferred name
  # used by beyondcell and the MoA).
  info <- FindDrugs(bc, x = rownames(bc@scaled), na.rm = FALSE)
  # Switch points.
  sp <- data.frame(switch.point = bc@switch.point[info$bc.names],
                   row.names = info$bc.names)
  # One plot per level.
  p4s <- lapply(lvl[in.lvl], function(l) {
    ### Subset residuals' means and switch points.
    res <- bc@ranks[[idents]][info$bc.names,
                              paste0("residuals.mean.", l), drop = FALSE]
    colnames(res) <- "residuals.mean"
    df <- transform(merge(res, sp, by = 0), row.names = Row.names, Row.names = NULL)
    ### Residual's deciles.
    res.decil <- quantile(as.numeric(res$residuals.mean),
                          prob = seq(from = 0, to = 1, length = 11))
    ### Drug annotation.
    sp_lower_01 <- as.numeric(df$switch.point) < 0.1
    sp_lower_06 <- as.numeric(df$switch.point) < 0.6
    sp_higher_04 <- as.numeric(df$switch.point) > 0.4
    sp_higher_09 <- as.numeric(df$switch.point) > 0.9
    res_lower_10 <- as.numeric(df$residuals.mean) < res.decil[["10%"]]
    res_higher_90 <- as.numeric(df$residuals.mean) > res.decil[["90%"]]
    df$annotation <- rep("no", times = nrow(df))
    df$annotation[sp_lower_01 & res_higher_90] <- "TOP-HighSensitivityDrugs"
    df$annotation[sp_higher_09 & res_lower_10] <- "TOP-LowSensitivityDrugs"
    df$annotation[sp_higher_04 & sp_lower_06 &
                    res_lower_10] <- "TOP-Differential-LowSensitivityDrugs"
    df$annotation[sp_higher_04 & sp_lower_06 &
                    res_higher_90] <- "TOP-Differential-HighSensitivityDrugs"
    ### Drug labels.
    df$labels <- rep(NA, times = nrow(df))
    decreasing_order <- c("TOP-Differential-HighSensitivityDrugs",
                          "TOP-HighSensitivityDrugs")
    unique.annotations <- unique(df$annotation[df$annotation != "no"])
    sel.labels <- unlist(sapply(unique.annotations, function(x) {
      sub.df <- subset(df, subset = df$annotation == x)
      if (x %in% decreasing_order) {
        sub.df <- sub.df[order(sub.df$residuals.mean, sub.df$switch.point,
                               decreasing = TRUE), ]
      } else {
        sub.df <- sub.df[order(sub.df$residuals.mean, sub.df$switch.point,
                               decreasing = FALSE), ]
      }
      return(rownames(sub.df)[1:min(top, nrow(sub.df))])
    }))
    df[sel.labels, "labels"] <- info$preferred.and.sigs[match(sel.labels,
                                                              info$bc.names)]
    ### Topnames.
    if(length(topnames[in.topnames]) > 0) {
      topnames <- FindDrugs(bc, x = topnames[in.topnames], na.rm = FALSE)
      df[match(topnames$bc.names,
               table = rownames(df)), "labels"] <- topnames$preferred.and.sigs
    }
    ### Colours and names.
    colors <- c("#1D61F2", "#DA0078", "orange", "#C7A2F5", "grey80", "black")
    names <- c("TOP-LowSensitivityDrugs", "TOP-HighSensitivityDrugs",
               "TOP-Differential-HighSensitivityDrugs",
               "TOP-Differential-LowSensitivityDrugs", "no", "black")
    ### Circle's borders colour.
    df$borders <- df$annotation
    df$borders[df$labels != ""] <- "black"
    ### Reorder df so labeled points are plotted on top.
    df <- rbind(subset(df, subset = df$borders != "black"),
                subset(df, subset = df$borders == "black"))
    ### Plot.
    p <- ggplot(df, aes(x = as.numeric(residuals.mean),
                        y = as.numeric(switch.point), color = borders,
                        fill = annotation)) +
      geom_point(shape = 21, alpha = alpha, size = pt.size) +
      scale_color_manual(values = setNames(colors, names)) +
      scale_fill_manual(values = setNames(colors, names), breaks = names[1:4],
                        drop = FALSE) + theme_classic() +
      geom_vline(xintercept = res.decil["10%"], linetype = "dotted") +
      geom_vline(xintercept = res.decil["90%"], linetype = "dotted") +
      geom_hline(yintercept = 0.9, linetype = "dotted") +
      geom_hline(yintercept = 0.1, linetype = "dotted") +
      geom_hline(yintercept = 0.4, linetype = "dotted") +
      geom_hline(yintercept = 0.6, linetype = "dotted") + ylim(0, 1) +
      labs(title = paste(idents, "=", lvl),
           caption = paste0("x cut-offs: first and last deciles; y cut-offs:",
                            " 0.1, 0.4, 0.6 and 0.9")) + 
      xlab("Residuals' Mean") + ylab("Switch Point") +
      ggrepel::geom_text_repel(label = df$labels, force = force, na.rm = TRUE,
                               ...) +
      guides(fill = guide_legend(title = "Drug Annotation"), color = FALSE) +
      cowplot::theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
    return(p)
  })
  return(p4s)
}
