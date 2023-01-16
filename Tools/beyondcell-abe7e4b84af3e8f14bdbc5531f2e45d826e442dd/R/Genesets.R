#' @title Creates a geneset object
#' @description This function creates a \code{\link[beyondcell]{geneset}}
#' object.
#' @name GenerateGenesets
#' @importFrom qusage read.gmt
#' @importFrom gdata trim
#' @param x A pre-loaded matrix, a ranked matrix or a path to a GMT file with
#' custom gene sets. See Details for more information.
#' @param n.genes Number of up and/or down-regulated genes used to compute each
#' signature.
#' @param mode Whether the output \code{geneset} must contain up and/or
#' down-regulated genes. See Details for more information.
#' @param filters If \code{x} is a pre-loaded matrix, you can provide a list of
#' filters to subset it. You can specify which \code{drugs}, sig \code{IDs},
#' mechanisms of action (\code{MoAs}), \code{targets} and/or \code{sources} you
#' are interested in (cap insensitive). You can call
#' \code{\link[beyondcell]{ListFilters}} to check all the available values for
#' these filters. The signatures that pass \strong{ANY} of them are included in
#' the output.
#' @param comparison \code{"treated_vs_control"} or
#' \code{"sensitive_vs_resistant"}. See Details for more information.
#' @param include.pathways Logical. Return \code{beyoncell}'s pre-computed
#' signatures for functional pathways?
#' @details \code{x} can be:
#' \itemize{
#' \item{A pre-loaded matrix:} {Either \code{PSc}, \code{SSc} or \code{DSS}.}
#' \item{A ranked matrix:} {A matrix with genes as rows and signatures as
#' columns that contains some type of numeric value such a t-stat or a LFC to
#' rank the genes accordingly.}
#' \item{A path to a GMT file:} {A file that contains custom gene sets. Each
#' gene set must have an "_UP" or "_DOWN" suffix.}
#' }
#' In addition, \code{mode} can be:
#' \itemize{
#' \item{\code{"up"}:} {To compute the signatures using only up-regulated
#' genes.}
#' \item{\code{"down"}:} {To compute the signatures using only down-regulated
#' genes.}
#' \item{\code{c("up", "down")} :} {To compute the signatures using both up and
#' down-regulated genes.}
#' }
#' If \code{x} is a path to a GMT file, \code{mode} is deprecated and the names
#' of all gene sets must end in "_UP" or "_DOWN" to indicate the \code{mode} of
#' each one.
#'
#' Finally, \code{comparison} can be:
#' \itemize{
#' \item{\code{"treated_vs_control"}:} {(\code{PSc} and \code{DSS} like) When
#' the numeric values or the gene sets in the GMT file were obtained from a
#' comparison between drug treated and untreated cells.}
#' \item{\code{"sensitive_vs_resistant"}:} {(\code{SSc} like) When the numeric
#' values or the gene sets in the GMT file were obtained from a comparison
#' between drug sensitive and resistant cells.}
#' }
#' When \code{x} is a pre-loaded matrix, \code{comparison} is set automatically.
#' @return A \code{geneset} object.
#' @examples
#' @export

GenerateGenesets <- function(x, n.genes = 250, mode = c("up", "down"),
                             filters = list(drugs = NULL, IDs = NULL, MoAs = NULL,
                                            targets = NULL, sources = NULL),
                             comparison = NULL, include.pathways = TRUE) {
  # --- Global Checks ---
  # Check if x is a pre-loaded matrix, a ranked matrix or a path to a GMT file.
  is.D <- c(identical(x, PSc), identical(x, SSc), identical(x, DSS))
  if (any(is.D)) {
    type <- "pre-loaded matrix"
    message(paste('Reading', c("PSc", "SSc", "DSS")[is.D], 'signatures...'))
    n.max <- 500
  } else if (is.matrix(x) | is.data.frame(x)) {
    type <- "matrix"
    message('Reading input matrix...')
    ### Check if x is numeric.
    if (!is.numeric(x)) stop('x must be a numeric matrix.')
    x <- as.matrix(x)
    ### Check if there are missing values.
    if (sum(is.na(x)) > 0) {
      warning('x contains NAs that will not be used to compute the geneset object.')
    }
    ### Check if there are duplicated values in the same column.
    dup.values <- apply(x, 2, FUN = function(y) any(duplicated(na.omit(y))))
    if (is.null(colnames(x))) colnames(x) <- 1:ncol(x)
    if(is.null(rownames(x))) rownames(x) <- 1:nrow(x)
    if (any(dup.values)) {
      warning(paste0('The following columns contain duplicated values:',
                     paste0(colnames(x)[dup.values], collapse = ", "), '.'))

    }
    n.max <- max(apply(x, 2, FUN = function(z) length(na.omit(z))))
  } else {
    type <- "gmt"
    message('Reading gmt file...')
    gmt.file <- tryCatch(qusage::read.gmt(x), error = function(cond) {
      message(paste0('The GMT file does not seem to exist: ', x, ' .'))
      stop(paste('x must be either a pre-loaded matrix (PSc, SSc or DSS), a
                 ranked matrix or the path to a GMT file.'))
    }, warning = function(cond) {
      message(paste('Warning when reading the GMT file. Here is the',
                    'original warning message:'))
      warning(paste0(cond))
      return(suppressWarnings(qusage::read.gmt(x)))
    })
    ### Check for duplicated gene sets.
    upper.gmt.names <- toupper(names(gmt.file))
    if (anyDuplicated(upper.gmt.names) != 0) {
      duplicated.gene.sets <- unique(names(gmt.file)[duplicated(upper.gmt.names)])
      stop(paste0('The GMT file contains duplicated gene set\'s: ',
                  paste0(duplicated.gene.sets, collapse = ", "), '.'))
    }
  }
  # Check n.genes and mode.
  if (any(!(mode %in% c("up", "down")))) stop('Incorrect mode.')
  mode <- sort(unique(mode), decreasing = TRUE)
  if (type != "gmt") {
    ### Number of genes.
    if (length(n.genes) != 1 | n.genes[1]%%1 != 0 | n.genes[1] < 1) {
      stop('n.genes must be a positive integer.')
    } else if (n.genes > n.max) {
      stop(paste0('n.genes exceeds the maximum number of genes in signature (',
                  n.max, ').'))
    }
  } else {
    ### Number of genes.
    if (!identical(n.genes, 250)) warning('x is a GMT file, n.genes is deprecated.')
    ### Mode in GMT files.
    n.up <- length(unique(grep(pattern = "_UP$", x = upper.gmt.names)))
    n.down <- length(unique(grep(pattern = "_DOWN$", x = upper.gmt.names)))
    if (n.up + n.down != length(names(gmt.file))) {
      stop('All gene sets\' names in the GMT file must end in "_UP" or "_DOWN".')
    } else {
      if (n.up > 0 & n.down > 0) {
        if (!identical(mode, c("up", "down"))) {
          mode <- c("up", "down")
          warning(paste('The GMT file includes UP and DOWN gene sets. mode',
                        'changed to c("up", "down").'))
        }
      } else if (n.up > 0) {
        if (mode != "up") {
          mode <- "up"
          warning('The GMT file only includes UP gene sets. mode changed to "up".')
        }
      } else if (n.down > 0) {
        if (mode != "down") {
          mode <- "down"
          warning('The GMT file only includes DOWN gene sets. mode changed to "down".')
        }
      }
    }
  }
  # Check filters and comparison.
  filters.names <- c("drugs", "IDs", "MoAs", "targets", "sources")
  selected.filters <- names(filters)
  if (any(!(selected.filters %in% filters.names))) stop('Invalid names in filters.')
  if (type != "pre-loaded matrix") {
    ### Filters.
    filters_class <- sapply(filters, is.null)
    if (any(!filters_class)) {
      warning('x is not a pre-loaded matrix, filters is deprecated.')
    }
    ### Comparison.
    if (is.null(comparison)) {
      stop(paste('Comparison must be either "treated_vs_control" or',
                 '"control_vs_treated".'))
    } else if (length(comparison) != 1 |
               !(comparison[1] %in% c("treated_vs_control", "control_vs_treated"))) {
      stop(paste('Comparison must be either "treated_vs_control" or',
                 '"control_vs_treated".'))
    }
  } else {
    ### Filters.
    filters_class <- sapply(filters, is.null) | sapply(filters, is.character)
    if (any(!filters_class)) {
      stop(paste0('Incorrect value for filter\'s entry: "',
                  paste0(selected.filters[!filters_class], collapse = ", "),
                  '". You must provide a character vector.'))
    }
    selected.filters <- selected.filters[!sapply(filters, is.null)]
    ### Comparison.
    if (is.null(comparison)) {
      if(is.D[2]) comparison <- "control_vs_treated"
      else comparison <- "treated_vs_control"
    } else {
      if (length(comparison) != 1 |
          !(comparison[1] %in% c("treated_vs_control", "control_vs_treated"))) {
        stop('Incorrect comparison.')
      }
      if (is.D[2] & comparison != "control_vs_treated") {
        comparison <- "control_vs_treated"
        warning('x = SSc, comparison changed to "control_vs_treated".')
      } else if (!is.D[2] & comparison != "treated_vs_control") {
        comparison <- "treated_vs_control"
        warning(paste0('x = ', c("PSc", "SSc", "DSS")[is.D], ', comparison ',
                       'changed to "treated_vs_control".'))
      }
    }
  }
  # Check include.pathways.
  if (length(include.pathways) != 1 | !is.logical(include.pathways)) {
    stop('include.pathways must be TRUE or FALSE.')
  }
  # --- Code ---
  # If x is a pre-loaded matrix...
  if (type == "pre-loaded matrix") {
    ### sig IDs.
    if (is.D[1]) {
      info <- subset(drugInfo, subset = drugInfo$sources == "LINCS")
    } else if (is.D[2]) {
      info <- subset(drugInfo, subset = drugInfo$sources != "LINCS")
    } else if (is.D[3]) {
      info <- subset(drugInfo, subset = drugInfo$IDs %in% DSS[[1]]$sig_id)
      x <- PSc # DSS is a subset of PSc
    }
    if (length(selected.filters) == 0) {
      ids <- unique(info$IDs)
    } else {
      ids <- unique(unlist(lapply(selected.filters, function(y) {
        tryCatch(suppressWarnings(GetIDs(values = filters[[y]], filter = y,
                                         df = info)),
                 error = function(cond) character())
      })))
      warnings <- unlist(lapply(selected.filters, function(z) {
        tryCatch(GetIDs(values = filters[[z]], filter = z, df = info),
                 error = function(cond) {
                   err <- paste0(z, ": ", paste0(filters[[z]],
                                                 collapse = ", "), ".\n")
                   return(err)
                 }, warning = function(cond) {
                   warn <- as.character(cond)
                   warn.values <- strsplit(sapply(strsplit(warn, split = ": "),
                                                  `[[`, 3), split = ", ")
                   return(paste0(z, ": ", warn.values))
                 })
      }))
      warnings <- warnings[!startsWith(warnings, prefix = "sig_")]
      if (length(ids) == 0) {
        stop('Couldn\'t find signatures that matched any of the filters.')
      } else if (length(warnings) > 0) {
        warning(paste('The following filters\' values yielded no results:\n',
                      paste0("   - ", warnings, " ", collapse = "")))
      }
    }
    genes <- lapply(ids, function(sig) {
      l <- list(up = x[["up"]][1:n.genes, sig], down = x[["down"]][1:n.genes, sig])
      return(l)
    })
    names(genes) <- ids
    # Else if x is a ranked matrix...
  } else if (type == "matrix") {
    genes <- apply(x, 2, FUN = function(sig) {
      l <- list()
      if ("up" %in% mode) {
        up <- na.omit(rownames(x)[order(sig, decreasing = TRUE,
                                        na.last = NA)[1:n.genes]])
        l <- c(l, list(up = up))
      }
      if ("down" %in% mode) {
        down <- na.omit(rownames(x)[order(sig, decreasing = FALSE,
                                          na.last = NA)[1:n.genes]])
        l <- c(l, list(down = down))
      }
      return(l)
    })
    # Else if x is a GMT file.
  } else if (type == "gmt") {
    unique.gene.sets <- unique(gsub(pattern = "_UP$|_DOWN$", replacement = "",
                                    x = names(gmt.file), ignore.case = TRUE))
    genes <- setNames(lapply(unique.gene.sets, function(sig) {
      l <- list()
      if (toupper(paste0(sig, "_UP")) %in% upper.gmt.names) {
        l <- c(l, list(up = gmt.file[[match(toupper(paste0(sig, "_UP")),
                                            table = upper.gmt.names)]]))
      }
      if (toupper(paste0(sig, "_DOWN")) %in% upper.gmt.names) {
        l <- c(l, list(down = gmt.file[[match(toupper(paste0(sig, "_DOWN")),
                                              table = upper.gmt.names)]]))
      }
      return(l)
    }), unique.gene.sets)
  }
  # Drug IDs.
  if (type == "pre-loaded matrix") {
    info <- subset(info, subset = info$IDs %in% ids)
    info <- aggregate(.~ IDs, data = info, na.action = NULL, FUN = function(rw) {
      paste(na.omit(unique(rw)), collapse = ", ")
    })
    info <- info[order(info$IDs, decreasing = FALSE), ]
  } else {
    info <- data.frame()
  }
  # Pathways.
  if (include.pathways) {
    paths <- lapply(pathways, function(p) p[names(p)[mode %in% names(p)]])
  } else {
    paths <- list()
  }
  # Output.
  return(geneset(genelist = c(genes, paths), n.genes = n.genes,
                 mode = mode, info = info, comparison = comparison))
}

#' @title Returns all the possible values for the specified filter
#' @description This function returns all the available values for \code{drugs},
#' \code{IDs}, \code{MoAs}, \code{targets} or \code{sources} filters in
#' \code{\link[beyondcell]{GenerateGenesets}} function.
#' @name ListFilters
#' @param entry Either \code{"drugs"}, \code{"IDs"}, \code{"MoAs"},
#' \code{"targets"} or \code{"sources"}.
#' @return All the possible values for the specified \code{entry}.
#' @export

ListFilters <- function(entry) {
  # --- Checks and Code ---
  if (entry == "drugs") {
    out <- sort(unique(drugInfo$drugs), decreasing = FALSE)
  } else if (entry == "IDs") {
    out <- sort(unique(drugInfo$IDs), decreasing = FALSE)
  } else if (entry == "MoAs") {
    out <- sort(unique(drugInfo$MoAs), decreasing = FALSE)
  } else if (entry == "targets") {
    out <- sort(unique(drugInfo$targets), decreasing = FALSE)
  } else if (entry == "sources") {
    out <- sort(unique(drugInfo$sources), decreasing = FALSE)
  } else {
    stop("Incorrect entry.")
  }
  return(out)
}
