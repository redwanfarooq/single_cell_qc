# QC helper functions

# Get 10x count matrix
get.10x.matrix <- function(file, type = NULL, remove.suffix = FALSE) {
  matrix <- Matrix::readMM(file) %>% as("CsparseMatrix")
  barcodes <- vroom::vroom(file.path(dirname(file), "barcodes.tsv.gz"), delim = "\t", col_names = FALSE)
  features <- vroom::vroom(file.path(dirname(file), "features.tsv.gz"), delim = "\t", col_names = FALSE)
  colnames(matrix) <- barcodes[, 1, drop = TRUE]
  rownames(matrix) <- features[, 1, drop = TRUE]
  if (remove.suffix) rownames(matrix) <- gsub(pattern = "-1$", replacement = "", rownames(matrix))
  if (!is.null(type)) matrix <- matrix[features[, 3, drop = TRUE] == type, ]
  return(matrix)
}

# Get barcounter count matrix
get.barcounter.matrix <- function(file, include.total = FALSE, features.pattern = NULL, add.suffix = FALSE) {
  matrix <- vroom::vroom(file, delim = ",", col_names = TRUE) %>%
    tibble::column_to_rownames("cell_barcode") %>%
    as.matrix() %>%
    t()
  if (!include.total) matrix <- matrix[-1, ]
  if (!is.null(features.pattern)) matrix <- matrix[grepl(pattern = features.pattern, rownames(matrix)), ]
  if (add.suffix) colnames(matrix) <- paste(colnames(matrix), "1", sep = "-")
  return(matrix)
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
  data <- data %>% as.data.frame() %>% 
    dplyr::group_by({{ pos }}, {{ fill }}) %>% 
    dplyr::summarise(count = dplyr::n())

  plot <- ggplot2::ggplot(data = data) +
    ggplot2::geom_col(mapping = ggplot2::aes(x = count, y = {{ pos }}, fill = {{ fill }}), position = position_stack(reverse = TRUE)) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal()
  return(plot)
}

# Generate HTO log-counts density ridge plots
logcounts.plot <- function(matrix, pseudocount = 1, samples = NULL, title = NULL) {
  data <- matrix %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(across(everything(), ~ log10(.x + pseudocount)))
  if (!is.null(samples)) data <- dplyr::bind_cols(data, sample = samples)
  data <- tidyr::pivot_longer(data, cols = starts_with("HTO"), names_to = "HTO", values_to = "count")

  plot <- ggplot2::ggplot(data = data)
  if (is.null(samples)) {
    plot <- plot +
      ggridges::geom_density_ridges(mapping = ggplot2::aes(x = count, y = HTO), fill = "#dddddd", show.legend = FALSE)
  } else {
    plot <- plot +
      ggridges::geom_density_ridges(
        mapping = ggplot2::aes(x = count, y = sample, fill = sample),
        show.legend = FALSE
      ) +
      ggplot2::facet_wrap(~HTO, scales = "free_x", ncol = 2)
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

# Add cell metadata to SingleCellExperiment
add.cell.metadata <- function(sce, metadata) {
  if (!identical(colnames(sce), rownames(metadata))) {
    warning("Cell barcodes do not match, retaining only barcodes in common")
    common <- intersect(colnames(sce), rownames(metadata))
    sce <- sce[, common]
    metadata <- metadata[common, , drop = FALSE]
  }
  df <- colData(sce) %>%
    merge(metadata, by = "row.names")
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  colData(sce) <- df %>% as("DataFrame")
  return(sce)
}

# Generate UMAP plot of known and predicted multiplets
multiplets.plot <- function(sce) {
  sce$Type <- dplyr::case_when(
    sce$known ~ "known",
    sce$predicted ~ "predicted",
    TRUE ~ "singlet"
  )

  plot <- scater::ggcells(sce, mapping = ggplot2::aes(x = UMAP.1, y = UMAP.2, colour = Type)) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::scale_colour_manual(values = c("singlet" = "grey", "known" = "#d82526", "predicted" = "#ffc156")) +
    ggplot2::labs(title = "Multiplets", x = "UMAP_1", y = "UMAP_2") +
    ggplot2::theme_minimal()
  return(plot)
}
