#' @title Computes a UMAP projection and therapeutic clusters using the BCS
#' @description This function uses the beyondcell scores (BCS) to compute a
#' UMAP projection of the data and to cluster cells according to their
#' sensitivity to the tested drugs (therapeutic clusters).
#' @name bcUMAP
#' @import scales
#' @import Seurat
#' @import ggplot2
#' @param bc \code{\link[beyondcell]{beyondcell}} object.
#' @param pc Number of principal components to use (which can be estimated from
#' an elbow plot). If \code{pc = NULL} (default), the function will stop prior
#' to compute the UMAP projection and the therapeutic clusters. See Details for
#' more information.
#' @param k.neighbors (\code{\link[Seurat]{FindNeighbors}}' \code{k.param})
#' Defines \emph{k} for the k-Nearest Neighbour algorithm.
#' @param res (\code{\link[Seurat]{FindClusters}}' \code{resolution}) Value of
#' the resolution parameter, use a value above/below 1.0 if you want to obtain
#' a larger/smaller number of communities. Can be a single number or a numeric
#' vector.
#' @param add.DSS Use background BCS computed with \code{DSS} signatures
#' (\code{add.DSS = TRUE}) or just use the signatures included in the \code{bc}
#' object (\code{add.DSS = FALSE}) to compute the UMAP projection and the
#' therapeutic clusters. If the number of drugs in \code{bc}
#' (excluding pathways) is <= 20, it is recomended to set \code{add.DSS = TRUE}.
#' Note that if \code{add.DSS = TRUE}, the regression and subset steps that have
#' been applied on \code{bc} will also be applied on the background BCS.
#' @param method (\code{\link[Seurat]{RunUMAP}}'s \code{umap.method}) UMAP
#' implementation to run. Can be:
#' \itemize{
#' \item{\code{uwot}}: Runs the UMAP via the \code{\link[uwot]{umap}} function
#' of the \code{uwot} package.
#' \item{\code{umap-learn}}: Runs the UMAP via the \code{Seurat} wrapper of the
#' python \code{umap-learn} package.
#' }
#' @param return.model (\code{\link[Seurat]{RunUMAP}}'s \code{return.model})
#' Whether \code{RunUMAP} will return the \code{uwot} model.
#' @param elbow.path Path to save the elbow plot. If \code{elbow.path = NULL}
#' (default), the plot will be printed.
#' @details This function performs all the steps required to obtain a UMAP
#' reduction of the data and cluster the cells according to the BCS.
#'
#' You will normally require to run the function twice:
#' \enumerate{
#' \item Using \code{pc = NULL} to obtain the elbow plot.
#' \item Specifying a value for the \code{pc} parameter according to this plot.
#' This second time, the UMAP reduction and the therapeutic clusters will be
#' computed.}
#'
#' Note that \code{add.DSS} must be the same in both runs, so the elbow plot
#' obtained in 1 is still valid in 2. If \code{add.DSS = TRUE}, the background
#' BCS will be stored in the \code{bc} object and the function will skip this
#' step the second time.
#' @return A \code{beyondcell} object with the UMAP reduction in
#' \code{@@reductions} slot and the therapeutic clusters for each \code{res}
#' in \code{bc@@meta.data}. Also, an elbow plot (\code{\link[ggplot2]{ggplot2}}
#' object) is printed (if \code{elbow.path = NULL}) or saved (if
#' \code{elbow.path != NULL}).
#' @examples
#' @export

bcUMAP <- function(bc, pc = NULL, k.neighbors = 20, res = 0.2,
                   add.DSS = FALSE, method = "uwot", return.model = FALSE,
                   elbow.path = NULL) {
  # --- Checks ---
  # Check that bc is a beyondcell object.
  if (class(bc) != "beyondcell") stop('bc must be a beyondcell object.')
  # Check pc.
  if (!is.null(pc)) {
    if (length(pc) != 1 | pc[1] < 2) stop('pc must be an integer >= 2.')
  }
  # Check k.neighbors.
  if (length(k.neighbors) != 1 | k.neighbors < 1) {
    stop('k.neighbors must be a positive integer.')
  }
  # Check res.
  if (any(sapply(res, function(x) !is.numeric(x))) |
      any(sapply(res, function(y) y < 0))) {
    stop('res must be a vector of numbers >= 0.')
  }
  # Check add.DSS.
  not.paths <- !(rownames(bc@normalized) %in% names(pathways))
  n.drugs <- sum(not.paths)
  if (length(add.DSS) != 1 | !is.logical(add.DSS)) {
    stop('add.DSS must be TRUE or FALSE.')
  } else if (!add.DSS) {
    if (n.drugs <= 10) {
      stop(paste('Only', n.drugs, 'drug signatures (excluding pathways) are',
                 'present in the bc object, please set add.DSS = TRUE.'))
    } else if (n.drugs <= 20) {
      warning(paste('Computing an UMAP reduction for', n.drugs,
                    'drugs. We recommend to set add.DSS = TRUE when the number',
                    'of signatures (excluding pathways) is below or equal to 20.'))
    }
  }
  # Check method.
  if (length(method) != 1 | !(method[1] %in% c("uwot", "umap-learn"))) {
    stop('Incorrect method.')
  }
  # Check return.model
  if (length(return.model) != 1 | !is.logical(return.model)) {
    stop('return.model must be TRUE or FALSE.')
  }
  if (method == "umap-learn" & return.model == TRUE) {
    warning(paste('return.model = TRUE is only valid when method =',
                  '"umap-learn". Changing return.model to FALSE.'))
    return.model <- FALSE
  }
  # Check elbow.path.
  if (!is.null(elbow.path)) {
    if (length(elbow.path) != 1 | !is.character(elbow.path)) {
      stop(paste('To save the elbow plot, you must specify a character string',
                 'with the desired destination.'))
    } else {
      parent.dir <- sub(pattern = "(.*\\/)([^.]+)(\\.[[:alnum:]]+$)",
                        replacement = "\\1", x = elbow.path)
      if (!dir.exists(parent.dir)) {
        stop(paste(parent.dir, 'path not found.'))
      }
    }
  }
  # --- Code ---
  # Cells in bc.
  cells <- colnames(bc@normalized)
  if (add.DSS) {
    ### DSS (background) BCS.
    if (!identical(sort(rownames(bc@background), decreasing = FALSE),
                   sort(DSS[[1]]$sig_id, decreasing = FALSE)) |
        !identical(sort(colnames(bc@background), decreasing = FALSE),
                   sort(cells, decreasing = FALSE)) |
        !identical(bc@regression$order, bc@regression$order.background)) {
      message('Computing background BCS using DSS signatures...')
      ### Genesets.
      gs.background <- suppressMessages(
        GenerateGenesets(DSS, n.genes = bc@n.genes, mode = bc@mode,
                         include.pathways = FALSE))
      ### BCS.
      background <- suppressWarnings(
        bcScore(bc@expr.matrix, gs = gs.background, expr.thres = bc@thres))
      ### Add metadata.
      background@meta.data <- background@meta.data[, -c(1:ncol(background@meta.data))]
      background <- bcAddMetadata(background, metadata = bc@meta.data)
      ### Subset and regress (if needed).
      if (bc@regression$order[1] == "subset") {
        background <- bcSubset(background, cells = cells)
      } else if (bc@regression$order[1] == "regression") {
        message('Regressing background BCS...')
        background <- suppressMessages(
          bcRegressOut(background, vars.to.regress = bc@regression[["vars"]]))
      }
      if (bc@regression$order[2] == "subset") {
        background <- bcSubset(background, cells = cells)
      } else if (bc@regression$order[2] == "regression") {
        message('Regressing background BCS...')
        background <- suppressMessages(
          bcRegressOut(background, vars.to.regress = bc@regression[["vars"]]))
      }
      ### Add background@normalized to bc@background.
      bc@background <- background@normalized
      ### Add order.background to bc@regression.
      bc@regression[["order.background"]] <- bc@regression[["order"]]
    } else {
      message('Background BCS already computed. Skipping this step.')
    }
    ### Add background to bc.
    all.rows <- unique(c(rownames(bc@normalized), rownames(bc@background)))
    merged.score <- rbind(bc@normalized, bc@background[, cells])[all.rows, ]
    ### Scale.
    merged.score <- t(apply(merged.score, 1, scales::rescale, to = c(0, 1)))
    bc.merged <- beyondcell(scaled = merged.score)
  } else {
    ### No background BCS.
    message(paste('DSS background not computed. UMAP will be created just with',
                  'the drugs (not pathways) in bc object.'))
    bc.merged <- bc
  }
  sc <- Seurat::CreateSeuratObject(bc.merged@scaled[not.paths, , drop = FALSE])
  # PCA.
  sc <- Seurat::ScaleData(sc, features = rownames(sc), do.scale = FALSE,
                          do.center = FALSE)
  sc <- Seurat::RunPCA(sc, features = rownames(sc), npcs = 100,  maxit = 100000)
  # Elbow plot.
  elbowplot <- Seurat::ElbowPlot(sc, ndims = 50) +
    ggplot2::theme(legend.position = "bottom")
  if (is.null(elbow.path)) {
    message('Printing elbow plot...')
    print(elbowplot)
  } else {
    message(paste('Saving elbow plot in', elbow.path))
    ggplot2::ggsave(elbow.path, plot = elbowplot)
  }
  if (!is.null(pc)) {
    ### Clusters.
    message("Obtaining therapeutic clusters...")
    sc <- Seurat::FindNeighbors(sc, dims = 1:pc, k.param = k.neighbors)
    sc <- Seurat::FindClusters(sc, resolution = res)
    ### UMAP.
    message('Computing beyondcell\'s UMAP reduction...')
    sc <- Seurat::RunUMAP(sc, dims = 1:pc, umap.method = method,
                          return.model = return.model, n.components = 2,
                          verbose = FALSE)
    bc@reductions <- sc@reductions
    ### Therapeutic clusters.
    message('Adding therapeutic clusters to metadata...')
    meta <- sc@meta.data[, paste0("RNA_snn_res.", res), drop = FALSE]
    colnames(meta) <- paste0("bc_clusters_res.", res)
    ### bc@metadata without repeated beyondcell clusters.
    repeated.res <- colnames(bc@meta.data) %in% colnames(meta)
    bc.metadata <- bc@meta.data[, !repeated.res, drop = FALSE]
    ### bc@meta.data rownames.
    meta.cells <- rownames(bc@meta.data)
    ### Merge bc@meta.data (without repeated beyondcell clusters) and meta. If
    ### the bc object was subsetted and some cells removed, their bc_clusters
    ### will be NA).
    new.meta <- transform(merge(bc.metadata, meta, by = 0, all.x = TRUE),
                          row.names = Row.names, Row.names = NULL)[meta.cells, ]
    bc@meta.data <- new.meta
  }
  return(bc)
}
