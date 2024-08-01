# ==============================
# Helper functions
# ==============================


# Explicitly define pipe operator
`%>%` <- magrittr::`%>%`


# C++ function for efficiently reading BarCounter CSV file to sparse matrix
Rcpp::sourceCpp(
  code = r"(
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
List read_barcounter_csv(string file) {
  // read header line and extract feature names
  ifstream infile(file);
  string line;
  getline(infile, line);
  istringstream iss(line);
  string feature;
  vector<string> features;
  while(getline(iss, feature, ',')) {
    features.emplace_back(feature);
  }
  features.erase(features.begin()); // remove barcode column header

  // initialise vectors to store barcodes, row indices, values, and column pointers
  vector<string> barcodes;
  vector<int> i, x;
  vector<int> p = {0};

  // read remaining lines and extract barcodes, row indices, values, and column pointers
  while (getline(infile, line)) {
    // extract barcode
    istringstream iss(line);
    string barcode;
    getline(iss, barcode, ',');
    barcodes.emplace_back(barcode);

    // read values
    vector<int> y;
    string val;
    while(getline(iss, val, ',')) {
      y.emplace_back(stoi(val));
    }

    // extract row indices and values for non-zero entries
    for (size_t j = 0; j < y.size(); ++j) {
      if (y[j] > 0) {
        i.emplace_back(j);
        x.emplace_back(y[j]);
      }
    }

    // set next column pointer
    p.emplace_back(i.size());
  }

  return List::create(
    Named("i") = i,
    Named("p") = p,
    Named("x") = x,
    Named("features") = features,
    Named("barcodes") = barcodes
  );
}
)"
)


#' Modified implementation of CellRanger OrdMag algorithm
#'
#' Finds the `q`th quantile of `x` and returns a threshold `r`-fold less
#' than this value.
#'
#' @param x Numeric vector of log10-transformed counts.
#' @param q Numeric scalar between 0 and 1 (default 0.99).
#' @param r Numeric scalar greater than 1 (default 10).
#'
#' @returns Numeric scalar.
.ordmag <- function(x, q = 0.99, r = 10) {
  return(quantile(x, q, na.rm = TRUE) - log10(r))
}


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
  if (!all(is.list(user), !is.null(names(user)))) stop("'user' must be a named list")
  if (!all(is.list(default), !is.null(names(default)))) stop("'default' must be a named list")

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


#' Prepend index to a string
#'
#' Utility function for generating dynamic numbered lists in Markdown. Searches for
#' a variable named `list.index` in the calling environment which is used to keep a
#' running count of the index (automatically incremented with each call).
#'
#' @param x Character scalar specifying list item.
#' @param i Integer scalar specifying index.
#'
#' @returns A character scalar of `x` prepended with `i`.
#'
#' @export
prepend.index <- function(x) {
  i <- get("list.index", envir = parent.frame()) %>% as.integer()
  str <- glue::glue("{i}. {x}")
  list.index <<- i + 1
  return(str)
}


#' Get 10x count matrix
#' 
#' Loads 10x count matrix (in HDF5 format).
#' 
#' @param file Path to file.
#' @param version Character scalar. 10x HDF5 version ('auto', 'v2', or 'v3').
#' @param cells Character vector. If specified, will subset to matching cell
#' barcodes.
#' @param type Character vector. If specified, will subset to features with matching
#' values in feature type field (v3 only).
#' @param remove.suffix Logical scalar (default `TRUE`). Remove '-1' suffix automatically
#' appended to cell barcodes by Cell Ranger.
#' 
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#' as column names.
#'
#' @export
get.10x.h5 <- function(file,
                       version = c("auto", "v2", "v3"),
                       cells = NULL,
                       type = NULL,
                       remove.suffix = TRUE,
                       group = NULL) {
  if (!file.exists(file)) stop(file, " does not exist")

  infile <- hdf5r::H5File$new(filename = file, mode = "r")

  version <- match.arg(version)
  if (version == "auto") {
    version <- if (hdf5r::existsGroup(infile, "matrix")) "v3" else "v2"
    message("Detected 10x HDF5 version: ", version)
  }
  features <- if (version == "v2") "genes" else "features/id"

  if (is.null(group)) {
    if (version == "v3") {
      group <- "matrix"
    } else {
      group <- names(infile)
      if (length(group) > 1) stop("Multiple groups detected: ", paste(group, collapse = ", "), ". Please specify 'group' parameter.")
    }
  }
  counts <- infile[[paste(group, "data", sep = "/")]]
  indices <- infile[[paste(group, "indices", sep = "/")]]
  indptr <- infile[[paste(group, "indptr", sep = "/")]]
  shape <- infile[[paste(group, "shape", sep = "/")]]
  features <- infile[[paste(group, features, sep = "/")]]
  barcodes <- infile[[paste(group, "barcodes", sep = "/")]]
  matrix <- Matrix::sparseMatrix(
    i = indices[],
    p = indptr[],
    x = as.numeric(counts[]),
    dims = shape[],
    dimnames = list(features[], barcodes[]),
    index1 = FALSE,
    repr = "C"
  )

  if (remove.suffix) colnames(matrix) <- gsub(pattern = "-1$", replacement = "", colnames(matrix))

  if (!is.null(cells)) matrix <- matrix[, match(cells, colnames(matrix))]
  if (any(is.na(colnames(matrix)))) {
    warning(sum(is.na(colnames(matrix))), " barcode(s) in 'cells' not present in count matrix")
    matrix <- matrix[, !is.na(colnames(matrix))]
  }

  if (version == "v3") {
    feature.type <- infile[["matrix/features/feature_type"]]
    if (!is.null(type)) {
      matrix <- matrix[feature.type[] %in% type, ]
    } else {
      if (length(unique(feature.type[])) > 1) message("Multiple feature types detected: ", paste(unique(feature.type[]), collapse = ", "), ". Returning list of matrices; please specify 'type' parameter to return a single matrix.")
      matrix <- lapply(
        unique(feature.type[]),
        function(type, matrix, feature.type) matrix[grep(pattern = type, x = feature.type), ],
        matrix = matrix,
        feature.type = feature.type[]
      ) |>
        setNames(unique(feature.type[]))
    }

    return(matrix)
  }
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
#' @param remove.suffix Logical scalar (default `TRUE`). Remove '-1' suffix automatically
#'  appended to cell barcodes by Cell Ranger.
#'
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#'  as column names.
#'
#' @export
get.10x.matrix <- function(file,
                           cells = NULL,
                           type = NULL,
                           remove.suffix = TRUE) {
  if (!file.exists(file)) stop(file, " does not exist")

  matrix <- Matrix::readMM(file) %>% as("CsparseMatrix")
  barcodes <- get.10x.barcodes(file.path(dirname(file), "barcodes.tsv.gz"), remove.suffix = remove.suffix)
  features <- get.10x.features(file.path(dirname(file), "features.tsv.gz"))
  colnames(matrix) <- rownames(barcodes)
  rownames(matrix) <- rownames(features)
  if (!is.null(cells)) matrix <- matrix[, match(cells, colnames(matrix))]
  if (any(is.na(colnames(matrix)))) {
    warning(sum(is.na(colnames(matrix))), " barcode(s) in 'cells' not present in count matrix")
    matrix <- matrix[, !is.na(colnames(matrix))]
  }
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
  if (!file.exists(file)) stop(file, " does not exist")

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
  if (!file.exists(file)) stop(file, " does not exist")

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
#' @param cells Character vector. If specified, will subset to matching cell
#'  barcodes (after adding suffix if `add.suffix = TRUE`).
#' @param features.pattern Regular expression string. If specified, with subset
#'  to features with matching names.
#' @param include.total Logical scalar (default `FALSE`). Include column containing
#'  total UMI count per barcode.
#' @param add.suffix Logical scalar (default `FALSE`). Add '-1' suffix to match
#'  cell barcodes from Cell Ranger.
#'
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#'  as column names.
#'
#' @export
get.barcounter.matrix <- function(file,
                                  cells = NULL,
                                  features.pattern = NULL,
                                  include.total = FALSE,
                                  add.suffix = FALSE) {
  if (!file.exists(file)) stop(file, " does not exist")

  # Efficiently read BarCounter CSV file into sparse matrix without loading dense matrix into memory
  sparse <- read_barcounter_csv(file)
  matrix <- Matrix::sparseMatrix(
    i = sparse$i,
    p = sparse$p,
    x = sparse$x,
    dims = c(length(sparse$features), length(sparse$barcodes)),
    dimnames = list(sparse$features, sparse$barcodes),
    index1 = FALSE,
    repr = "C"
  )
  if (add.suffix) colnames(matrix) <- paste(colnames(matrix), "1", sep = "-")
  if (!include.total) matrix <- matrix[rownames(matrix) != "total", ]
  if (!is.null(cells)) matrix <- matrix[, match(cells, colnames(matrix))]
  if (any(is.na(colnames(matrix)))) {
    warning(sum(is.na(colnames(matrix))), " barcode(s) in 'cells' not present in count matrix")
    matrix <- matrix[, !is.na(colnames(matrix))]
  }
  if (!is.null(features.pattern)) matrix <- matrix[grepl(pattern = features.pattern, x = rownames(matrix)), ]

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
add.cell.metadata <- function(x,
                              metadata,
                              replace = c("matching", "all", "none")) {
  if (!any(is(x, "SingleCellExperiment"), is(x, "Seurat"))) {
    stop("x must be an object of class 'SingleCellExperiment' or 'Seurat'")
  }
  replace <- match.arg(replace)

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
add.feature.metadata <- function(x,
                                 metadata,
                                 assay = NULL,
                                 replace = c("matching", "all", "none")) {
  if (!any(is(x, "SingleCellExperiment"), is(x, "Seurat"))) {
    stop("x must be an object of class 'SingleCellExperiment' or 'Seurat'")
  }
  replace <- match.arg(replace)

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


#' Modified implementation of CellRanger empirical expected cell number estimation
#' algorithm
#'
#' Estimates the expected cell number *x* from empirical total UMI counts per barcode
#' by minimising the loss function (OrdMag(*x*) - *x*)^2 / *x* over the range of *x*
#' between 10 and 45000 in steps of 10, where OrdMag(*x*) is the number of barcodes
#' with UMI count above a threshold 10-fold less than the 99th percentile of the top
#' *x* barcodes.
#'
#' @param counts Integer vector of total UMI counts per barcode.
#'
#' @returns Integer scalar *x*.
#'
#' @export
estimate.expected.cells <- function(counts) {
  df <- data.frame(
    x = seq_along(counts),
    logcounts = sort(log10(counts), decreasing = TRUE),
    stat = NA
  ) %>%
    filter(x <= 45000)

  for (i in seq.int(from = 10, to = nrow(df), by = 10)) {
    threshold <- .ordmag(df$logcounts[seq_len(i)])
    df$stat[i] <- (sum(df$logcounts > threshold) - df$x[i])^2 / df$x[i]
  }

  return(df$x[match(min(df$stat, na.rm = TRUE), df$stat)])
}


#' Multimodal cell calling algorithm
#'
#' @description
#' Modified implementation of CellRanger ARC cell calling algorithm.
#'
#' Algorithm steps:
#' 1. Compute log10-transformed total counts for each modality per barcode
#' 2. Filter to barcodes with count >= 1 for every modality
#' 3. Deduplicate barcodes with identical counts across all modalities
#' 4. Determine initial thresholds for cell-containing barcodes using OrdMag algorithm
#'    for each modality and compute centroids for cell-containing and empty barcodes
#' 5. Peform k-means clustering (k = 2) initialised using centroids from OrdMag-based
#'    thresholds and update barcode labels based on cluster assignment
#' 6. Project barcode labels from deduplicated set to full set
#'
#' Key differences from CellRanger ARC implementation are:
#' * Lack of initial pre-filtering steps
#' * Generalisation to more than 2 modalities
#'
#' @param matrix.list List of raw count matrices.
#' @param n.expected.cells Integer scalar. Targeted cell recovery.
#' @param ordmag.quantile Numeric scalar (default 0.99). Quantile used for OrdMag function.
#' @param ordmag.ratio Numeric scalar (default 10). Ratio used for OrdMag function.
#'
#' @returns Character vector of cell-containing barcodes.
#'
#' @export
multimodal.cell.caller <- function(matrix.list,
                                   n.expected.cells,
                                   ordmag.quantile = 0.99,
                                   ordmag.ratio = 10) {
  if (!is.list(matrix.list)) stop("'matrix.list' must be a list")
  if (length(matrix.list) < 2) stop("'matrix.list' must contain at least 2 matrices")
  if (is.null(names(matrix.list))) names(matrix.list) <- as.character(seq_along(matrix.list))

  matrix.list <- lapply(matrix.list, function(x) x[, colSums(x) > 0])

  logcounts <- Map(
    f = function(x, i) {
      data.frame(log10(colSums(x))) %>%
        setNames(paste("logcounts", i, sep = ".")) %>%
        tibble::rownames_to_column("barcode")
    },
    x = matrix.list,
    i = names(matrix.list)
  ) %>%
    Reduce(
      x = .,
      f = function(x, y) dplyr::full_join(x, y, by = "barcode")
    ) %>%
    dplyr::filter(dplyr::if_all(dplyr::everything(), function(x) !is.na(x)))

  deduplicated <- logcounts %>%
    dplyr::select(dplyr::starts_with("logcounts")) %>%
    dplyr::distinct()
  centers <- deduplicated %>%
    dplyr::mutate(
      cell = dplyr::if_else(
        dplyr::if_all(
          dplyr::everything(),
          function(x) {
            threshold <- .ordmag(
              x = sort(x, decreasing = TRUE)[seq_len(n.expected.cells)],
              q = ordmag.quantile,
              r = ordmag.ratio
            )
            return(x > threshold)
          }
        ),
        TRUE,
        FALSE
      )
    ) %>%
    dplyr::group_by(cell) %>%
    dplyr::summarise(dplyr::across(dplyr::starts_with("logcounts"), mean)) %>%
    dplyr::arrange(cell) %>%
    tibble::column_to_rownames("cell") %>%
    as.matrix()
  deduplicated <- deduplicated %>%
    dplyr::mutate(
      cell = dplyr::case_match(kmeans(x = ., centers = centers)$cluster, 1 ~ FALSE, 2 ~ TRUE)
    )

  cells <- dplyr::left_join(
    x = logcounts,
    y = deduplicated,
    by = paste("logcounts", names(matrix.list), sep = ".")
  ) %>%
    dplyr::filter(cell) %>%
    dplyr::pull("barcode")

  return(cells)
}


#' Combine results of multiple QC filters
#'
#' Takes a list or data frame of QC filtering results (logical vectors with
#' 1 element per cell barcode) and combines the results to identify barcodes which
#' have passed all filters.
#'
#' @param filters List or data frame of logical vectors.
#' @param missing Logical scalar (default FALSE). Value to convert NA values
#' in `filters` to prior to combining results.
#'
#' @returns Logical vector.
#'
#' @export
combine.filters <- function(filters, missing = FALSE) {
  if ((!is.list(filters) || !is.data.frame(filters)) && !all(sapply(filters, is.logical))) stop("filters must be a list of logical vectors")

  for (i in seq_along(filters)) {
    filters[[i]] <- ifelse(is.na(filters[[i]]), missing, filters[[i]])
  }

  return(Reduce(magrittr::or, filters, init = FALSE))
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
#' @param title Character scalar specifying plot title.
#' @param x.lab Character scalar specifying x-axis label.
#' @param y.lab Character scalar specifying y-axis label.
#' @param colour.lab Character scalar specifying colour legend label.
#' @param log Log-transform axis scales:
#'  * `TRUE` to transform all axes
#'  * Character vector with values 'x' and/or 'y' to transform specified axes
#'
#' @returns Plot object.
#'
#' @export
scatter.plot <- function(data,
                         x,
                         y,
                         colour,
                         ...,
                         title = NULL,
                         x.lab = NULL,
                         y.lab = NULL,
                         colour.lab = NULL,
                         log = NULL) {
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
violin.plot <- function(data,
                        x,
                        y,
                        colour,
                        ...,
                        points = TRUE,
                        title = NULL,
                        x.lab = NULL,
                        y.lab = NULL,
                        colour.lab = NULL,
                        log = NULL) {
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
#' @param title Character scalar specifying plot title.
#' @param x.lab Character scalar specifying x-axis label.
#'
#' @returns Plot object.
#'
#' @export
hist.plot <- function(data,
                      x,
                      ...,
                      title = NULL,
                      x.lab = NULL) {
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
#' @param title Character scalar specifying plot title.
#' @param pos.lab Character scalar specifying categorical axis label.
#' @param fill.lab Character scalar specifying fill legend label.
#' @param log Log-transform axis scales:
#'  * `TRUE` to transform all axes
#'  * Character vector with values 'x' and/or 'y' to transform specified axes
#'
#' @returns Plot object.
#'
#' @export
bar.plot <- function(data,
                     pos,
                     fill = pos,
                     ...,
                     title = NULL,
                     pos.lab = NULL,
                     fill.lab = NULL,
                     log = NULL) {
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
#' as cells using the output of [DropletUtils::barcodeRanks()].
#'
#' @param bcrank.res DataFrame output of [DropletUtils::barcodeRanks()].
#' @param cells Character vector of cell-containing barcodes.
#' @param ... Arguments to pass to [scatter.plot()]
#'
#' @returns Plot object.
#'
#' @export
bcrank.plot <- function(bcrank.res, cells, ...) {
  data <- bcrank.res %>%
    as.data.frame() %>%
    dplyr::mutate(cell = dplyr::if_else(rownames(.) %in% cells, TRUE, FALSE)) %>%
    dplyr::group_by(rank, total, cell) %>%
    dplyr::summarise(fraction.cells = sum(cell) / dplyr::n())

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


#' Generate feature log-counts distribution plots
#'
#' Makes log-counts distribution plots for selected features, optionally
#' grouping by sample.
#'
#' @param matrix Matrix or matrix-like object containing counts with
#'  features as row names.
#' @param ... Fixed aesthetics to pass to [ggridges::geom_density_ridges()].
#' @param features.pattern Character vector of regular expressions to select
#'  features for plotting.
#' @param pseudocount Numeric scalar (default 1). Added to counts before log
#'  transformation (to avoid infinite values).
#' @param x.lower Numeric scalar (default 0.3). Lower limit of x-axis scale
#'  (choose a value greater than 0 to avoid loss of resolution due to inflated
#'  counts close to zero; set to `NA` to use ggplot2 default).
#' @param x.upper Numeric scalar (default maximum log-count in `matrix`).
#'  Upper limit of x-axis scale (set to `NA` for to use ggplot2 default).
#' @param samples Character vector of the same length as number of
#'  cells specifying sample of origin. If specified, will group by sample.
#' @param title Character scalar specifying plot title.
#'
#' @returns Plot object.
#'
#' @export
logcounts.plot <- function(matrix,
                           ...,
                           features.pattern = "^[a-zA-Z]+_",
                           pseudocount = 1,
                           x.lower = 0.3,
                           x.upper = max(log10(matrix + pseudocount)),
                           samples = NULL,
                           title = NULL) {
  data <- matrix %>%
    as.matrix() %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(across(everything(), ~ log10(.x + pseudocount)))
  if (!is.null(samples)) data <- dplyr::bind_cols(data, Sample = samples)
  data <- tidyr::pivot_longer(data, cols = matches(features.pattern), names_to = "Feature", values_to = "count")

  plot <- ggplot2::ggplot(data = data)
  if (is.null(samples)) {
    plot <- plot + ggridges::geom_density_ridges(mapping = ggplot2::aes(x = count, y = Feature), ...)
  } else {
    plot <- plot +
      ggridges::geom_density_ridges(mapping = ggplot2::aes(x = count, y = Sample, fill = Sample), ...) +
      ggplot2::facet_wrap(~Feature, scales = "free_x", ncol = 2, labeller = "label_both")
  }
  plot <- plot +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::xlim(x.lower, x.upper) +
    ggplot2::labs(
      title = title,
      x = sprintf("log10(UMI count + %d)", pseudocount),
      y = ifelse(is.null(samples), "Feature", "Sample")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0))

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
multiplet.plot <- function(x,
                           rd.res,
                           ...,
                           projection = c("TSNE", "UMAP"),
                           use.dimred = NULL,
                           random.seed = NULL) {
  if (!is(x, "SingleCellExperiment")) stop("x must be an object of class 'SingleCellExperiment'")

  x <- add.cell.metadata(x, rd.res)
  x$Type <- dplyr::case_when(
    x$known ~ "known",
    x$predicted ~ "predicted",
    TRUE ~ "singlet"
  )

  projection <- match.arg(projection)
  if (!projection %in% names(SingleCellExperiment::reducedDims(x))) {
    set.seed(random.seed)
    if (projection == "TSNE") {
      x <- scater::runTSNE(x, dimred = use.dimred)
    } else if (projection == "UMAP") {
      x <- scater::runUMAP(x, dimred = use.dimred)
    }
  }

  var.x <- rlang::sym(paste(projection, "1", sep = "."))
  var.y <- rlang::sym(paste(projection, "2", sep = "."))
  plot <- scater::ggcells(x, mapping = ggplot2::aes(x = {{ var.x }}, y = {{ var.y }}, colour = Type)) +
    ggplot2::geom_point(...) +
    ggplot2::scale_colour_manual(values = c("singlet" = "grey80", "known" = "#d82526", "predicted" = "#ffc156")) +
    ggplot2::labs(title = "Multiplets", x = paste(projection, "1", sep = "_"), y = paste(projection, "2", sep = "_")) +
    ggplot2::theme_minimal()

  return(plot)
}

#' Generate tSNE/UMAP plot of batch corrected cells
#'
#' Makes tSNE/UMAP plot of batch corrected cells using the output of
#' [batchelor::correctExperiments()].
#'
#' @inheritParams multiplet.plot
#' @param x SingleCellExperiment object output of
#'  [batchelor::correctExperiments()].
#' @param control Character scalar specifying name of control sample (i.e. a
#'  cross-batch technical replicate).
#' @param assay.type Character or integer scalar specifying which assay in `x`
#'  contains log-normalised expression values; only used if existing dimensionality
#'  reduction results not present in `reducedDims(x)`.
#' @param n.cells Integer scalar (default 10000). Specifies number of cells to
#'  plot; if `x` contains more cells, will randomly downsample to `n.cells`.
#'  Useful for reducing overplotting and reducing size of interactive plots. Set
#'  to `Inf` to force plotting all cells in `x`.
#' @param title Character scalar specifying plot title.
#'
#' @returns Plot object.
#'
#' @export
batch.plot <- function(x,
                       ...,
                       control = NULL,
                       assay.type = c("merged", "corrected"),
                       projection = c("TSNE", "UMAP"),
                       use.dimred = c("PCA", "corrected"),
                       n.cells = 10000,
                       random.seed = NULL,
                       title = NULL) {
  if (!is(x, "SingleCellExperiment")) stop("x must be an object of class 'SingleCellExperiment'")

  if (is.null(control)) {
    x$control <- FALSE
  } else {
    x$control <- ifelse(x$Sample == control, TRUE, FALSE)
  }

  if (length(SingleCellExperiment::reducedDims(x)) == 0) {
    set.seed(random.seed)
    if (!is.integer(assay.type)) assay.type <- match.arg(SummarizedExperiment::assayNames(x)[1], choices = assay.type)
    dec <- scran::modelGeneVar(x, assay.type = "logcounts", block = x$batch)
    hvg <- scran::getTopHVGs(dec, n = 5000)
    x <- scater::runPCA(x, assay.type = assay.type, subset_row = hvg)
  }

  if (!is.integer(use.dimred)) use.dimred <- match.arg(SingleCellExperiment::reducedDimNames(x), choices = use.dimred)
  projection <- match.arg(projection)
  if (!projection %in% names(SingleCellExperiment::reducedDims(x))) {
    set.seed(random.seed)
    if (projection == "TSNE") {
      x <- scater::runTSNE(x, dimred = use.dimred)
    } else if (projection == "UMAP") {
      x <- scater::runUMAP(x, dimred = use.dimred)
    }
  }

  var.x <- rlang::sym(paste(projection, "1", sep = "."))
  var.y <- rlang::sym(paste(projection, "2", sep = "."))
  plot <- scater::ggcells(
    x[, sample(seq_len(ncol(x)), min(n.cells, ncol(x)))],
    mapping = ggplot2::aes(x = {{ var.x }}, y = {{ var.y }}, fill = batch, colour = control)
  ) +
    ggplot2::geom_point(...) +
    ggplot2::labs(
      title = title,
      x = paste(projection, "1", sep = "_"),
      y = paste(projection, "2", sep = "_"),
      fill = "Batch",
      colour = "Control"
    ) +
    ggplot2::theme_minimal()

  return(plot)
}
