# QC helper functions

# Update R Markdown params list
update.params <- function(default, user) {
  if (is(default, "list")) {
    new <- vector("list", length(default)) %>% setNames(names(default))
    for (x in names(default)) new[[x]] <- update.params(default[[x]], user[[x]])
  } else {
    if (!is.null(user)) new <- user else new <- default
  }
  return(new)
}

# Get 10x count matrix
get.10x.matrix <- function(file, type = NULL, remove.suffix = FALSE) {
  matrix <- Matrix::readMM(file) %>% as("CsparseMatrix")
  barcodes <- get.10x.barcodes(file, remove.suffix = remove.suffix)
  features <- get.10x.features(file)
  colnames(matrix) <- rownames(barcodes)
  rownames(matrix) <- rownames(features)
  if (!is.null(type)) matrix <- matrix[features$Type == type, ]
  return(matrix)
}

# Get 10x barcode metadata
get.10x.barcodes <- function(file, remove.suffix = FALSE) {
  df <- vroom::vroom(
    file.path(dirname(file), "barcodes.tsv.gz"),
    delim = "\t",
    col_names = "Barcode",
    show_col_types = FALSE
  ) %>%
    as.data.frame()
  if (remove.suffix) df$Barcode <- gsub(pattern = "-1$", replacement = "", df$Barcode)
  rownames(df) <- df$Barcode
  return(df)
}

# Get 10x features metadata
get.10x.features <- function(file, type = NULL) {
  df <- vroom::vroom(
    file.path(dirname(file), "features.tsv.gz"),
    delim = "\t",
    col_names = c("ID", "Symbol", "Type", "Chr", "Start", "End"),
    show_col_types = FALSE
  ) %>%
    as.data.frame()
  rownames(df) <- df$ID
  if (!is.null(type)) df <- df[df$Type == type, ]
  return(df)
}

# Get barcounter count matrix
get.barcounter.matrix <- function(file, include.total = FALSE, features.pattern = NULL, add.suffix = FALSE) {
  matrix <- vroom::vroom(file, delim = ",", col_names = TRUE, show_col_types = FALSE) %>%
    tibble::column_to_rownames("cell_barcode") %>%
    as.matrix() %>%
    t()
  if (!include.total) matrix <- matrix[-1, ]
  if (!is.null(features.pattern)) matrix <- matrix[grepl(pattern = features.pattern, rownames(matrix)), ]
  if (add.suffix) colnames(matrix) <- paste(colnames(matrix), "1", sep = "-")
  return(matrix)
}

# Add cell metadata to SingleCellExperiment or Seurat object
add.cell.metadata <- function(x, metadata, replace = FALSE) {
  if (all(!is(x, "SingleCellExperiment"), !is(x, "Seurat"))) stop("x must be an object of class 'SingleCellExperiment' or 'Seurat'")
  if (!identical(colnames(x), rownames(metadata))) {
    warning("Cell barcodes do not match, retaining only barcodes in common")
    common <- intersect(colnames(x), rownames(metadata))
    x <- x[, common]
    metadata <- metadata[common, , drop = FALSE]
  }
  if (replace) {
    df <- metadata
  } else {
    if (is(x, "SingleCellExperiment")) {
      current <- colData(x)
    } else if (is(x, "Seurat")) {
      current <- x@meta.data
    }
    df <- current %>%
      merge(metadata, by = "row.names")
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
  }
  if ((is(x, "SingleCellExperiment"))) {
    colData(x) <- df %>% as("DataFrame")
  } else if (is(x, "Seurat")) {
    x@meta.data <- df %>% as("data.frame")
  }
  return(x)
}

# Add feature metadata to SingleCellExperiment
add.feature.metadata <- function(x, metadata, replace = FALSE) {
  if (!identical(rownames(x), rownames(metadata))) {
    warning("Features do not match, retaining only features in common")
    common <- intersect(rownames(x), rownames(metadata))
    x <- x[common, ]
    metadata <- metadata[common, ]
  }
  if (replace) {
    df <- metadata
  } else {
    df <- rowData(x) %>%
      merge(metadata, by = "row.names")
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
  }
  rowData(x) <- df %>% as("DataFrame")
  return(x)
}

# Generate barcode rank plot
bcrank.plot <- function(bcrank.res, ed.res) {
  data <- bcrank.res %>%
    as.data.frame() %>%
    cbind(FDR = ed.res$FDR) %>%
    dplyr::group_by(rank, total) %>%
    dplyr::summarise(fraction.cells = sum(FDR <= 0.001, na.rm = TRUE) / dplyr::n())

  plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = rank, y = total, colour = fraction.cells)) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_gradient2(low = "#d82526", mid = "#ffc156", high = "#69b764", midpoint = 0.5) +
    ggplot2::labs(title = "Barcode rank plot", x = "Rank", y = "Total UMI count", colour = "Fraction cells") +
    ggplot2::theme_minimal()
  return(plot)
}

# Generate Monte Carlo p-value histogram for background barcodes
pvalue.plot <- function(ed.res) {
  data <- ed.res %>%
    as.data.frame() %>%
    dplyr::filter(Total <= metadata(ed.res)$lower & Total > 0)

  plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_histogram(mapping = ggplot2::aes(x = PValue), binwidth = 0.01) +
    ggplot2::labs(title = "Monte Carlo p-value distribution for background barcodes", x = "p-value") +
    ggplot2::theme_minimal()
  return(plot)
}

# Generate bar plot
bar.plot <- function(data, pos, fill = pos, title = NULL) {
  data <- data %>%
    as.data.frame() %>%
    dplyr::group_by({{ pos }}, {{ fill }}) %>%
    dplyr::summarise(count = dplyr::n())

  plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_col(mapping = ggplot2::aes(x = count, y = {{ pos }}, fill = {{ fill }}), position = position_stack(reverse = TRUE)) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal()
  return(plot)
}

# Generate feature log-counts density ridge plots
logcounts.plot <- function(matrix, feature.pattern = "^HTO_", pseudocount = 1, samples = NULL, title = NULL) {
  data <- matrix %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(across(everything(), ~ log10(.x + pseudocount)))
  if (!is.null(samples)) data <- dplyr::bind_cols(data, sample = samples)
  data <- tidyr::pivot_longer(data, cols = matches(feature.pattern), names_to = "Feature", values_to = "count")

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

# Generate UMAP plot of known and predicted multiplets
multiplets.plot <- function(x) {
  x$Type <- dplyr::case_when(
    x$known ~ "known",
    x$predicted ~ "predicted",
    TRUE ~ "singlet"
  )

  plot <- scater::ggcells(x, mapping = ggplot2::aes(x = UMAP.1, y = UMAP.2, colour = Type)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::scale_colour_manual(values = c("singlet" = "grey80", "known" = "#d82526", "predicted" = "#ffc156")) +
    ggplot2::labs(title = "Multiplets", x = "UMAP_1", y = "UMAP_2") +
    ggplot2::theme_minimal()
  return(plot)
}

# Generate violin + beeswarm plot of QC metrics
violin.plot <- function(data, x, y, colour, title = NULL, x.lab = NULL, y.lab = NULL, log = NULL) {
  data <- data %>%
    as.data.frame()

  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }})) +
    ggplot2::geom_violin(colour = "grey50") +
    ggbeeswarm::geom_quasirandom(mapping = ggplot2::aes(colour = {{ colour }}), size = 0.5) +
    scale_colour_manual(values = c("grey80", "#d82526"), guide = "none") +
    labs(title = title, x = x.lab, y = y.lab) +
    theme_minimal()
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- "y"
    if ("x" %in% log) stop("Unable to log transform categorical x-axis scale")
    if ("y" %in% log) plot <- plot + scale_y_log10()
  }
  return(plot)
}

# Generate scatter plot of QC metrics
scatter.plot <- function(data, x, y, colour, title = NULL, x.lab = NULL, y.lab = NULL, log = NULL) {
  data <- data %>%
    as.data.frame()

  plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_point(mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, colour = {{ colour }}), size = 0.5) +
    scale_colour_manual(values = c("grey80", "#d82526"), guide = "none") +
    labs(title = title, x = x.lab, y = y.lab) +
    theme_minimal()
  if (!is.null(log)) {
    if (is.logical(log) && log) log <- c("x", "y")
    if ("x" %in% log) plot <- plot + scale_x_log10()
    if ("y" %in% log) plot <- plot + scale_y_log10()
  }
  return(plot)
}
