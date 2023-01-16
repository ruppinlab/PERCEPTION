#' @title Geneset class
#' @description An object to represent signatures.
#' @slot genelist A list of drug signatures and functional \code{pathways} (if
#' specified) with up and/or down-regulated genes.
#' @slot n.genes Argument passed to \code{\link[beyondcell]{GenerateGenesets}}.
#' Number of up and/or down-regulated genes per signature.
#' @slot mode Argument passed to \code{GenerateGenesets}. Whether the
#' \code{geneset} contains up and/or down-regulated genes.
#' @slot info Dataframe with drug signatures information, including sig IDs,
#' drug names, MoAs, target genes and data sources (LINCS, CTRP, GDSC or CCLE).
#' This slot is only filled if \code{GenerateGenesets}' input is a pre-loaded
#' matrix.
#' @slot comparison Argument passed to \code{GenerateGenesets}. Either
#' \code{"treated_vs_control"} or \code{"control_vs_treated"}.

geneset <- setClass("geneset", slots = list(genelist = "list",
                                            n.genes = "numeric",
                                            mode = "character",
                                            info = "data.frame",
                                            comparison = "character"))

#' @title Beyondcell class
#' @description An object to represent the beyondcell scores (BCS) for each cell and
#' signature.
#' @slot scaled (Subsetted and/or regressed) scaled BCS.
#' @slot normalized (Subsetted and/or regressed) normalized BCS.
#' @slot data Original normalized BCS, without subsetting or regression.
#' @slot switch.point (Subsetted and/or regressed) scaled BCS for which the
#' normalized score in \code{@@data} is 0 (one switch point per signature).
#' @slot ranks List of dataframes with the BCS' statistics and ranks returned
#' by \code{\link[beyondcell]{bcRanks}}.
#' @slot expr.matrix Single-cell expression matrix used to compute the BCS.
#' @slot meta.data Dataframe that contains information about each cell
#' (including the therapeutic clusters and \code{\link[Seurat]{Seurat}}'s
#' \code{@@meta.data}).
#' @slot SeuratInfo List with information about the input \code{Seurat} object,
#' including the \code{@@reductions}.
#' @slot background (Subsetted and/or regressed) normalized BCS obtained using
#' DSS signatures. Useful to compute \code{beyondcell}'s UMAP reduction and the
#' therapeutic clusters when the number of drug signatures is low.
#' @slot reductions A list of dimensional reductions for this object.
#' @slot regression A list with the order of subset and regression steps
#' performed on the \code{beyondcell} object and the variables used for
#' regression.
#' @slot n.genes Argument passed to \code{\link[beyondcell]{GenerateGenesets}}.
#' Number of up and/or down-regulated genes per signature.
#' @slot mode Argument passed to \code{GenerateGenesets}. Whether the
#' \code{geneset} contains up and/or down-regulated genes.
#' @slot thres Argument \code{expr.thres} passed to
#' \code{\link[beyondcell]{bcScore}}. Minimum fraction of signature genes that
#' must be expressed in a cell to compute its BCS.

beyondcell <- setClass("beyondcell",
                       slots = list(scaled = "matrix", normalized = "matrix",
                                    data = "matrix", switch.point = "numeric",
                                    ranks = "list", expr.matrix = "matrix",
                                    meta.data = "data.frame", SeuratInfo = "list",
                                    background = "matrix", reductions = "list",
                                    regression = "list", n.genes = "numeric",
                                    mode = "character", thres = "numeric"))
