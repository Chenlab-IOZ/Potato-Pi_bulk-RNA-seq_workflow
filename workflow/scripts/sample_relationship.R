#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 15L) {
  stop(
    "Expected 15 arguments: counts, samples, three tables, three plots, ",
    "min_total_count, top_variable_genes, correlation_method, clustering_method, method, normalized_expression, normalization_factors"
  )
}

counts_path <- args[[1]]
samples_path <- args[[2]]
transformed_path <- args[[3]]
correlation_path <- args[[4]]
pca_scores_path <- args[[5]]
correlation_heatmap_path <- args[[6]]
pca_path <- args[[7]]
clustering_path <- args[[8]]
min_total_count <- as.numeric(args[[9]])
top_variable_genes <- as.integer(args[[10]])
correlation_method <- args[[11]]
clustering_method <- args[[12]]
method <- args[[13]]
normalized_path <- args[[14]]
factors_path <- args[[15]]
if (!method %in% c("deseq2", "edgeR")) {
  stop("method must be deseq2 or edgeR")
}
for (path in c(args[3:8], normalized_path, factors_path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

if (!is.finite(min_total_count) || min_total_count < 0) {
  stop("min_total_count must be a non-negative number")
}
if (is.na(top_variable_genes) || top_variable_genes < 2L) {
  stop("top_variable_genes must be an integer of at least 2")
}
if (!correlation_method %in% c("pearson", "spearman")) {
  stop("correlation_method must be pearson or spearman")
}
if (!clustering_method %in% c(
  "ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
  "median", "centroid"
)) {
  stop("Unsupported clustering_method: ", clustering_method)
}

counts_table <- read.delim(
  counts_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_metadata <- c("sample", "condition", "replicate", "batch")
metadata <- read.delim(
  samples_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(counts_table) < 3L || names(counts_table)[[1]] != "gene_id") {
  stop("Count matrix must start with gene_id and contain at least two samples")
}
if (anyDuplicated(counts_table$gene_id)) {
  stop("Count matrix contains duplicate gene_id values")
}
if (!all(required_metadata %in% names(metadata))) {
  stop("Sample metadata must contain: ", paste(required_metadata, collapse = ", "))
}
if (anyDuplicated(metadata$sample)) {
  stop("Sample metadata contains duplicate sample values")
}

sample_names <- names(counts_table)[-1]
if (!setequal(sample_names, metadata$sample)) {
  missing_metadata <- setdiff(sample_names, metadata$sample)
  missing_counts <- setdiff(metadata$sample, sample_names)
  stop(
    "Count matrix and metadata samples differ. Missing from metadata: ",
    paste(missing_metadata, collapse = ", "),
    "; missing from counts: ", paste(missing_counts, collapse = ", ")
  )
}
metadata <- metadata[match(sample_names, metadata$sample), required_metadata]
rownames(metadata) <- metadata$sample

count_matrix <- as.matrix(counts_table[, -1, drop = FALSE])
suppressWarnings(storage.mode(count_matrix) <- "numeric")
rownames(count_matrix) <- counts_table$gene_id
if (any(!is.finite(count_matrix))) {
  stop("Count matrix contains missing or non-numeric values")
}
if (any(count_matrix < 0) || any(count_matrix != round(count_matrix))) {
  stop("Count matrix must contain non-negative integer counts")
}
storage.mode(count_matrix) <- "integer"

keep <- rowSums(count_matrix) > min_total_count
count_matrix <- count_matrix[keep, , drop = FALSE]
if (nrow(count_matrix) < 2L) {
  stop("Fewer than two genes remain after count filtering")
}

if (any(colSums(count_matrix) == 0)) {
  stop("At least one sample has zero library size after filtering")
}
if (method == "edgeR") {
  dge <- edgeR::DGEList(counts = count_matrix)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  normalized_matrix <- edgeR::cpm(dge, log = FALSE)
  transformed_matrix <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  factors <- data.frame(
    sample = colnames(count_matrix),
    library_size = dge$samples$lib.size,
    TMM_factor = dge$samples$norm.factors,
    effective_library_size = dge$samples$lib.size * dge$samples$norm.factors
  )
  transformation <- "TMM log2-CPM (prior.count=2)"
} else {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix, colData = metadata, design = ~ 1
  )
  dds <- DESeq2::estimateSizeFactors(dds, type = "ratio")
  normalized_matrix <- DESeq2::counts(dds, normalized = TRUE)
  transformed_matrix <- SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
  factors <- data.frame(
    sample = colnames(count_matrix),
    library_size = colSums(count_matrix),
    size_factor = DESeq2::sizeFactors(dds)
  )
  transformation <- "DESeq2 rlog (blind=TRUE)"
}
write.table(
  data.frame(gene_id = rownames(normalized_matrix), normalized_matrix, check.names = FALSE),
  normalized_path, sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(factors, factors_path, sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  data.frame(
    gene_id = rownames(transformed_matrix),
    transformed_matrix,
    check.names = FALSE
  ),
  transformed_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_mad <- apply(transformed_matrix, 1, mad)
variable_genes <- names(sort(gene_mad[gene_mad > 0], decreasing = TRUE))
variable_genes <- head(variable_genes, top_variable_genes)
if (length(variable_genes) < 2L) {
  stop("Fewer than two variable genes remain for sample correlation")
}
correlation <- cor(
  transformed_matrix[variable_genes, , drop = FALSE],
  method = correlation_method
)
write.table(
  data.frame(sample = rownames(correlation), correlation, check.names = FALSE),
  correlation_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (any(!is.finite(correlation))) {
  stop("Non-finite sample correlations; check sample variation in selected genes")
}

annotation <- metadata[, c("condition", "batch"), drop = FALSE]
rownames(annotation) <- metadata$sample
correlation_distance <- as.dist(1 - correlation)
correlation_hclust <- hclust(correlation_distance, method = clustering_method)
pheatmap(
  correlation,
  cluster_rows = correlation_hclust,
  cluster_cols = correlation_hclust,
  annotation_row = annotation,
  annotation_col = annotation,
  display_numbers = TRUE,
  number_format = "%.2f",
  border_color = NA,
  main = paste0(
    method, " sample correlation (top ", length(variable_genes), " variable genes)"
  ),
  filename = correlation_heatmap_path,
  width = 12,
  height = 9
)

pca <- prcomp(t(transformed_matrix), center = TRUE, scale. = FALSE)
variance_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_scores <- data.frame(
  sample = rownames(pca$x),
  pca$x,
  metadata[rownames(pca$x), c("condition", "replicate", "batch")],
  check.names = FALSE,
  row.names = NULL
)
write.table(
  pca_scores,
  pca_scores_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pca_plot <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = condition, shape = batch, label = sample)
) +
  geom_point(size = 3.5) +
  geom_text_repel(
    show.legend = FALSE,
    seed = 17,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  labs(
    x = sprintf("PC1 (%.1f%%)", variance_percent[[1]]),
    y = sprintf("PC2 (%.1f%%)", variance_percent[[2]]),
    title = paste("PCA of samples:", transformation),
    color = "Condition",
    shape = "Batch"
  ) +
  theme_bw(base_size = 12) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(12, 20, 12, 20))
ggsave(pca_path, pca_plot, width = 9, height = 7)

sample_distance <- dist(t(transformed_matrix), method = "euclidean")
sample_hclust <- hclust(sample_distance, method = clustering_method)
condition_levels <- unique(metadata$condition)
condition_colors <- setNames(
  grDevices::hcl.colors(length(condition_levels), "Dark 3"),
  condition_levels
)
dendrogram <- as.dendrogram(sample_hclust)
dendrogram <- dendrapply(dendrogram, function(node) {
  if (is.leaf(node)) {
    sample <- attr(node, "label")
    node_par <- attr(node, "nodePar")
    node_par$lab.col <- condition_colors[[metadata[sample, "condition"]]]
    node_par$lab.cex <- 0.9
    attr(node, "nodePar") <- node_par
  }
  node
})

pdf(clustering_path, width = 10, height = 7)
par(mar = c(6, 4, 4, 2) + 0.1)
plot(dendrogram, main = paste("Hierarchical clustering:", method), ylab = "Euclidean distance")
legend(
  "topright",
  legend = condition_levels,
  col = condition_colors,
  pch = 15,
  title = "Condition",
  bty = "n"
)
dev.off()

cat("Input genes:", nrow(counts_table), "\n")
cat("Genes retained:", nrow(transformed_matrix), "\n")
cat("Variable genes used for correlation:", length(variable_genes), "\n")
cat("Samples:", ncol(transformed_matrix), "\n")

cat("Method:", method, "\n")
cat("Transformation:", transformation, "\n")
sessionInfo()
