# ==============================
# QC helper functions
# ==============================


#' Update R Markdown params list
#'
#' Recursively searches nested named list of default parameters for matching entries
#' in a user-defined list. If a match is found, the default value is replaced with
#' the user-defined value.
#'
#' @param user Named list containing user-defined parameters.
#' @param default Named list containing default parameters (named `params` by default
#'  in R Markdown documents).
#'
#' @returns A nested named list with the same structure as `default`.
#'
#' @export
update.params <- function(user, default = params) {
  new <- vector("list", length(default)) %>% setNames(names(default))
  for (x in names(default)) {
    if (x %in% names(user)) {
      if (is.list(default[[x]])) {
        new[[x]] <- update.params(user[[x]], default[[x]])
      } else {
        new[[x]] <- user[[x]]
      }
    } else {
      new[[x]] <- default[[x]]
    }
  }
  return(new)
}


#' Get 10x count matrix
#'
#' Loads 10x count matrix (in matrix market format).
#'
#' @param file Path to file.
#' @param cells Character vector. If specified, will subset to matching cell
#'  barcodes.
#' @param type Character vector. If specified, will subset to features with matching
#'  values in column 3 of 10x features TSV file.
#' @param remove.suffix Logical scalar (default `FALSE`). Remove '-1' suffix automatically
#'  appended to cell barcodes by Cell Ranger.
#'
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#'  as column names.
#'
#' @export
get.10x.matrix <- function(file, cells = NULL, type = NULL, remove.suffix = FALSE) {
  matrix <- Matrix::readMM(file) %>% as("CsparseMatrix")
  barcodes <- get.10x.barcodes(file.path(dirname(file), "barcodes.tsv.gz"), remove.suffix = remove.suffix)
  features <- get.10x.features(file.path(dirname(file), "features.tsv.gz"))
  colnames(matrix) <- rownames(barcodes)
  rownames(matrix) <- rownames(features)
  if (!is.null(cells)) matrix <- matrix[, cells]
  if (!is.null(type)) matrix <- matrix[features$Type %in% type, ]
  return(matrix)
}


#' Get 10x cell barcodes
#'
#' Loads 10x cell barcodes from TSV.
#'
#' @inheritParams get.10x.matrix
#'
#' @returns A data frame with the following columns:
#'  * Barcode: cell barcode
#'
#' @export
get.10x.barcodes <- function(file, remove.suffix = FALSE) {
  df <- vroom::vroom(
    file,
    delim = "\t",
    col_names = "Barcode",
    show_col_types = FALSE
  ) %>%
    as.data.frame()
  if (remove.suffix) df$Barcode <- gsub(pattern = "-1$", replacement = "", df$Barcode)
  rownames(df) <- df$Barcode
  return(df)
}


#' Get 10x feature metadata
#'
#' Loads 10x feature metadata from TSV.
#'
#' @inheritParams get.10x.matrix
#'
#' @returns A data frame with the following columns:
#'  * ID: feature ID
#'  * Symbol: feature symbol
#'  * Type: feature type
#'  * Chr: chromosome
#'  * Start: genomic start position
#'  * End: genomic end position
#'
#' @export
get.10x.features <- function(file, type = NULL) {
  df <- vroom::vroom(
    file,
    delim = "\t",
    col_names = c("ID", "Symbol", "Type", "Chr", "Start", "End"),
    show_col_types = FALSE
  ) %>%
    as.data.frame()
  rownames(df) <- df$ID
  if (!is.null(type)) df <- df[df$Type %in% type, ]
  return(df)
}


#' Get BarCounter count matrix
#'
#' Loads BarCounter count matrix (in CSV format).
#'
#' @param file Path to file.
#' @param include.total Logical scalar (default `FALSE`). Include column containing
#'  total UMI count per barcode.
#' @param cells Character vector. If specified, will subset to matching cell
#'  barcodes (after adding suffix if `add.suffix = TRUE`).
#' @param features.pattern Regular expression string. If specified, with subset
#'  to features with matching names.
#' @param add.suffix Logical scalar (default `FALSE`). Add '-1' suffix to match
#'  cell barcodes from Cell Ranger.
#'
#' @returns A matrix of counts with features as row names and cell barcodes
#'  as column names.
#'
#' @export
get.barcounter.matrix <- function(file, include.total = FALSE, cells = NULL, features.pattern = NULL, add.suffix = FALSE) {
  matrix <- vroom::vroom(file, delim = ",", col_names = TRUE, show_col_types = FALSE) %>%
    tibble::column_to_rownames("cell_barcode") %>%
    as.matrix() %>%
    t()
  if (add.suffix) colnames(matrix) <- paste(colnames(matrix), "1", sep = "-")
  if (!include.total) matrix <- matrix[-1, ]
  if (!is.null(cells)) matrix <- matrix[, cells]
  if (!is.null(features.pattern)) matrix <- matrix[grepl(pattern = features.pattern, rownames(matrix)), ]
  return(matrix)
}


#' Add cell metadata to SingleCellExperiment or Seurat object
#'
#' Merge or replace cell metadata contained in a SingleCellExperiment or
#' Seurat object with values in a data frame or data frame-like object.
#'
#' @param x A SingleCellExperiment or Seurat object.
#' @param metadata A data frame or data frame-like object.
#' @param replace Character scalar:
#'  * 'matching' (default): Replace metadata columns in `x` with
#'    matching columns in `metadata`.
#'  * 'all': Replace all metadata columns in `x` with columns in
#'    `metadata`.
#'  * 'none': Do not replace any metadata columns in `x`.
#'
#' @returns An object of the same class as `x` with updated cell metadata.
#'
#' @export
add.cell.metadata <- function(x, metadata, replace = c("matching", "all", "none")) {
  replace <- match.arg(replace)

  if (!any(is(x, "SingleCellExperiment"), is(x, "Seurat"))) stop("x must be an object of class 'SingleCellExperiment' or 'Seurat'")

  if (is(x, "SingleCellExperiment")) {
    current <- SummarizedExperiment::colData(x)
  } else if (is(x, "Seurat")) {
    current <- x@meta.data
  }

  if (!setequal(rownames(current), rownames(metadata))) stop("x and metadata must contain the same cell barcodes")
  metadata <- metadata[rownames(current), , drop = FALSE] # ensure matching row order

  if (replace == "matching") {
    current <- current[, !colnames(current) %in% intersect(colnames(current), colnames(metadata))]
  } else if (replace == "all") {
    current <- data.frame(row.names = rownames(current))
  }

  new <- cbind(current, metadata)

  if ((is(x, "SingleCellExperiment"))) {
    SummarizedExperiment::colData(x) <- new %>% as("DataFrame")
  } else if (is(x, "Seurat")) {
    x@meta.data <- new %>% as("data.frame")
  }
  return(x)
}


#' Add feature metadata to SingleCellExperiment or Seurat object
#'
#' Merge or replace feature metadata contained in a SingleCellExperiment or
#' Seurat object with values in a data frame or data frame-like object.
#'
#' @inheritParams add.cell.metadata
#' @param assay Name of assay containing relevant features (if `x` is a Seurat
#'  object). If not provided, will use `DefaultAssay(x)`.
#'
#' @returns An object of the same class as `x` with updated feature metadata.
#'
#' @export
add.feature.metadata <- function(x, metadata, assay = NULL, replace = c("matching", "all", "none")) {
  replace <- match.arg(replace)

  if (!any(is(x, "SingleCellExperiment"), is(x, "Seurat"))) stop("x must be an object of class 'SingleCellExperiment' or 'Seurat'")

  if (is(x, "SingleCellExperiment")) {
    current <- SummarizedExperiment::rowData(x)
  } else if (is(x, "Seurat")) {
    assay <- ifelse(is.null(assay), SeuratObject::DefaultAssay(x), assay)
    current <- x[[assay]]@meta.features
  }

  if (!setequal(rownames(current), rownames(metadata))) stop("x and metadata must contain the same features")
  metadata <- metadata[rownames(current), , drop = FALSE] # ensure matching row order

  if (replace == "matching") {
    current <- current[, !colnames(current) %in% intersect(colnames(current), colnames(metadata))]
  } else if (replace == "all") {
    current <- data.frame(row.names = rownames(current))
  }

  new <- cbind(current, metadata)

  if ((is(x, "SingleCellExperiment"))) {
    SummarizedExperiment::rowData(x) <- new %>% as("DataFrame")
  } else if (is(x, "Seurat")) {
    x[[assay]]@meta.features <- new %>% as("data.frame")
  }
  return(x)
}


#' Generate scatter plot
#'
#' Makes a scatter plot.
#'
#' @param data Data frame or object coercible to data frame.
#' @param x Column in `data` to determine x position.
#' @param y Column in `data` to determine y position.
#' @param colour Column in `data` to determine colour.
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_point()].
#' @param title Plot title.
#' @param x.lab x-axis label.
#' @param y.lab y-axis label.
#' @param colour.lab Colour legend label.
#' @param log Log-transform axis scales:
#'  * `TRUE` to transform all axes
#'  * Character vector with values 'x' and/or 'y' to transform specified axes
#'
#' @returns Plot object.
#'
#' @export
scatter.plot <- function(data, x, y, colour, ..., title = NULL, x.lab = NULL, y.lab = NULL, colour.lab = NULL, log = NULL) {
  data <- data %>%
    as.data.frame()

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, colour = {{ colour }})) +
    ggplot2::geom_point(...) +
    ggplot2::labs(title = title, x = x.lab, y = y.lab, colour = colour.lab) +
    ggplot2::theme_minimal()
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- c("x", "y")
    if ("x" %in% log) plot <- plot + ggplot2::scale_x_log10()
    if ("y" %in% log) plot <- plot + ggplot2::scale_y_log10()
  }
  return(plot)
}


#' Generate violin/beeswarm plot
#'
#' Makes a combined violin plot and beeswarm plot.
#'
#' @inheritParams scatter.plot
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_quasirandom()].
#' @param points Logical scalar (default `TRUE`). Add beeswarm points.
#'
#' @returns Plot object.
#'
#' @export
violin.plot <- function(data, x, y, colour, ..., points = TRUE, title = NULL, x.lab = NULL, y.lab = NULL, colour.lab = NULL, log = NULL) {
  data <- data %>%
    as.data.frame()

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, colour = {{ colour }})) +
    ggplot2::geom_violin(colour = "grey50") +
    ggplot2::labs(title = title, x = x.lab, y = y.lab, colour = colour.lab) +
    ggplot2::theme_minimal()
  if (points) {
    plot <- plot +
      ggbeeswarm::geom_quasirandom(...)
  }
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- "y"
    if ("x" %in% log) stop("Unable to log transform categorical x-axis scale")
    if ("y" %in% log) plot <- plot + ggplot2::scale_y_log10()
  }
  return(plot)
}


#' Generate histogram plot
#'
#' Makes a histogram plot.
#'
#' @param data Data frame or object coercible to data frame.
#' @param x Column in `data` to determine x position.
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_histogram()].
#' @param title Plot title.
#' @param x.lab x-axis label.
#'
#' @returns Plot object.
#'
#' @export
hist.plot <- function(data, x, ..., title = NULL, x.lab = NULL) {
  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }})) +
    ggplot2::geom_histogram(...) +
    ggplot2::labs(title = title, x = x.lab) +
    ggplot2::theme_minimal()
  return(plot)
}


#' Generate bar plot
#'
#' Makes a bar plot.
#'
#' @param data Data frame or object coercible to data frame.
#' @param pos Column in `data` to determine bar position.
#' @param fill Column in `data` to determine bar fill.
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_col()].
#' @param title Plot title.
#' @param pos.lab Categorical axis label.
#' @param fill.lab Fill legend label.
#' @param log Log-transform axis scales:
#'  * `TRUE` to transform all axes
#'  * Character vector with values 'x' and/or 'y' to transform specified axes
#'
#' @returns Plot object.
#'
#' @export
bar.plot <- function(data, pos, fill = pos, ..., title = NULL, pos.lab = NULL, fill.lab = NULL, log = NULL) {
  data <- data %>%
    as.data.frame() %>%
    dplyr::group_by({{ pos }}, {{ fill }}) %>%
    dplyr::summarise(count = dplyr::n())

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = count, y = {{ pos }}, fill = {{ fill }})) +
    ggplot2::geom_col(...) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(title = title, y = pos.lab, fill = fill.lab) +
    ggplot2::theme_minimal()
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- "x"
    if ("x" %in% log) plot <- plot + ggplot2::scale_x_log10()
    if ("y" %in% log) stop("Unable to log transform categorical y-axis scale")
  }
  return(plot)
}


#' Generate barcode rank plot
#'
#' Makes a barcode rank plot coloured by fraction of barcodes at each rank called
#' as cells using the output of [DropletUtils::barcodeRanks()] and
#' [DropletUtils::emptyDrops()].
#'
#' @param bcrank.res DataFrame output of [DropletUtils::barcodeRanks()].
#' @param ed.res DataFrame output of [DropletUtils::emptyDrops()].
#' @param ... Arguments to pass to [scatter.plot()]
#'
#' @returns Plot object.
#'
#' @export
bcrank.plot <- function(bcrank.res, ed.res, ...) {
  data <- bcrank.res %>%
    as.data.frame() %>%
    cbind(FDR = ed.res$FDR) %>%
    dplyr::group_by(rank, total) %>%
    dplyr::summarise(fraction.cells = sum(FDR <= 0.001, na.rm = TRUE) / dplyr::n())

  plot <- scatter.plot(
    data,
    x = rank,
    y = total,
    colour = fraction.cells,
    ...
  )
  return(plot)
}


#' Generate Monte Carlo p-value histogram for background barcodes
#'
#' Makes a histogram of computed p-values for droplets with UMI count less than
#' `lower` parameter used for [DropletUtils::emptyDrops()] (requires parameter
#' `test.ambient = TRUE`).
#'
#' @param ed.res DataFrame output of [DropletUtils::emptyDrops()].
#' @param ... Arguments to pass to [hist.plot()].
#'
#' @returns Plot object.
#'
#' @export
pvalue.plot <- function(ed.res, ...) {
  data <- ed.res %>%
    as.data.frame() %>%
    dplyr::filter(Total <= metadata(ed.res)$lower & Total > 0)

  plot <- hist.plot(
    data,
    x = PValue,
    ...
  )
  return(plot)
}


#' Generate feature log-counts density ridge plots
#'
#' Makes log-count density ridge plots for selected features, optionally
#' grouping by sample.
#'
#' @param matrix Matrix or matrix-like object containing counts with
#'  features as row names.
#' @param features.pattern Regular expression string to select features for
#'  plotting.
#' @param pseudocount Integer (default 1). Added to counts before log
#'  transformation (to avoid infinite values).
#' @param samples Character vector of the same length as number of
#'  cells specifying sample of origin. If specified, will group by sample.
#' @param title Plot title.
#'
#' @returns Plot object.
#'
#' @export
logcounts.plot <- function(matrix, features.pattern = "^[a-zA-Z]+_", pseudocount = 1, samples = NULL, title = NULL) {
  data <- matrix %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(across(everything(), ~ log10(.x + pseudocount)))
  if (!is.null(samples)) data <- dplyr::bind_cols(data, sample = samples)
  data <- tidyr::pivot_longer(data, cols = matches(features.pattern), names_to = "Feature", values_to = "count")

  plot <- ggplot2::ggplot(data = data)
  if (is.null(samples)) {
    plot <- plot +
      ggridges::geom_density_ridges(mapping = ggplot2::aes(x = count, y = Feature), fill = "grey80", show.legend = FALSE)
  } else {
    plot <- plot +
      ggridges::geom_density_ridges(
        mapping = ggplot2::aes(x = count, y = sample, fill = sample),
        show.legend = FALSE
      ) +
      ggplot2::facet_wrap(~Feature, scales = "free_x", ncol = 2)
  }
  plot <- plot +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::xlim(1, NA) +
    ggplot2::labs(title = title, x = sprintf("log10(UMI count + %d)", pseudocount), y = "density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(hjust = 0)
    )
  return(plot)
}


#' Generate tSNE/UMAP plot of known and predicted multiplets
#'
#' Makes tSNE/UMAP plot of known and predicted multiplets using the output of
#' [scDblFinder::recoverDoublets()].
#'
#' @param x A SingleCellExperiment object.
#' @param rd.res DataFrame output of [scDblFinder::recoverDoublets()].
#' @param ... Fixed aesthetics to pass to [ggplot2::geom_point()].
#' @param projection Character scalar. Determines type of projection used for
#'  plotting cells (will be computed if not present in `reducedDims(x)`):
#'  * 'TSNE': Use t-SNE projection
#'  * 'UMAP': Use UMAP projection
#' @param use.dimred Character or integer scalar specifying existing dimensionality
#'  reduction results to use for computing projection; only used if `projection`
#'  not present in `reducedDims(x)`.
#' @param random.seed Optional random seed for reproducibility of computed
#'  projection.
#'
#' @returns Plot object.
#'
#' @export
multiplet.plot <- function(x, rd.res, ..., projection = c("TSNE", "UMAP"), use.dimred = NULL, random.seed = NULL) {
  if (!is(x, "SingleCellExperiment")) stop("x must be an object of class 'SingleCellExperiment'")

  x <- add.cell.metadata(x, rd.res)

  projection <- match.arg(projection)
  var.x <- rlang::sym(paste(projection, "1", sep = "."))
  var.y <- rlang::sym(paste(projection, "2", sep = "."))

  x$Type <- dplyr::case_when(
    x$known ~ "known",
    x$predicted ~ "predicted",
    TRUE ~ "singlet"
  )

  if (!projection %in% names(SingleCellExperiment::reducedDims(x))) {
    set.seed(random.seed)
    if (projection == "TSNE") {
      x <- scater::runTSNE(x, dimred = use.dimred)
    } else if (projection == "UMAP") {
      x <- scater::runUMAP(x, dimred = use.dimred)
    }
  }

  plot <- scater::ggcells(x, mapping = ggplot2::aes(x = {{ var.x }}, y = {{ var.y }}, colour = Type)) +
    ggplot2::geom_point(...) +
    ggplot2::scale_colour_manual(values = c("singlet" = "grey80", "known" = "#d82526", "predicted" = "#ffc156")) +
    ggplot2::labs(title = "Multiplets", x = paste(projection, "1", sep = "_"), y = paste(projection, "2", sep = "_")) +
    ggplot2::theme_minimal()
  return(plot)
}
