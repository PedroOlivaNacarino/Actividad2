
# ============================================================
# RNA-seq – DESeq2 + Visualización FINAL
# ============================================================

# -------------------------------
# Librerías
# -------------------------------
library(tximport)
library(DESeq2)
library(readr)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(pheatmap)
library(reshape2)

# -------------------------------
# Importación de datos
# -------------------------------
files <- list.files("Quant", pattern="quant.sf", full.names=TRUE, recursive=TRUE)
names(files) <- gsub(".*/(.*)/quant.sf","\\1", files)

tx2gene <- read.delim("Transcrito_a_Gen.tsv", header=TRUE)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

design <- read.csv("Design.csv", row.names=1)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = design,
  design = ~ Condition
)

dds <- dds[rowSums(counts(dds)) > 1, ]

# -------------------------------
# DESeq2
# -------------------------------
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), "DESeq2_results.csv")

# -------------------------------
# Preparar resultados
# -------------------------------
res_df <- as.data.frame(res)
res_df <- rownames_to_column(res_df, var="Gene")

res_df <- res_df %>%
  mutate(Significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1,
                              "Yes", "No"))

sig_genes <- res_df %>% filter(Significant == "Yes")

# ============================================================
# VISUALIZACIÓN
# ============================================================


# ============================================================
# Volcano 
# ============================================================

volplot <- ggplot(res_df,
                  aes(x = log2FoldChange,
                      y = -log10(padj))) +
  
  # Puntos no significativos
  geom_point(data = subset(res_df, Significant == "No"),
             color = "grey75",
             size = 2,
             alpha = 0.6) +
  
  # Puntos significativos
  geom_point(data = subset(res_df, Significant == "Yes"),
             color = "#D55E00",
             size = 2.5,
             alpha = 0.9) +
  
  # Etiquetas de genes significativos
  geom_text_repel(
    data = sig_genes,
    aes(label = Gene),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  
  # Líneas de umbral
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "grey40",
             linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey40",
             linewidth = 0.5) +
  
  # Etiquetas y título
  labs(
    title = "Volcano plot – Expresión diferencial",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~p[adj])
  ) +
  
  # Tema visual
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "gray50", fill = NA, linewidth = 1),
    plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  )

volplot


# ============================================================
# MA plot 
# ============================================================

# Preparar datos (usar el res_df que ya tienes)
ma_df <- res_df %>%
  filter(!is.na(baseMean))

sig_ma <- ma_df %>% filter(Significant == "Yes")

MAplot <- ggplot(ma_df,
                 aes(x = log10(baseMean),
                     y = log2FoldChange)) +
  
  # Genes no significativos
  geom_point(data = subset(ma_df, Significant == "No"),
             color = "grey75",
             size = 2,
             alpha = 0.6) +
  
  # Genes significativos
  geom_point(data = subset(ma_df, Significant == "Yes"),
             color = "coral",
             size = 2.5,
             alpha = 0.9) +
  
  # Etiquetas de genes significativos
  geom_text_repel(
    data = sig_ma,
    aes(label = Gene),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  
  # Línea horizontal en 0
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey40",
             linewidth = 0.6) +
  
  # Límites similares al plotMA clásico
  coord_cartesian(ylim = c(-5, 5)) +
  
  # Etiquetas
  labs(
    title = "MA plot – Expresión diferencial",
    x = expression(log[10]~media~de~expresión),
    y = expression(log[2]~Fold~Change)
  ) +
  
  # Tema visual
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "gray50", fill = NA, linewidth = 1),
    plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.margin = unit(c(1,1,1,1), "lines")
  )

MAplot


# -------------------------------
# Heatmap
# -------------------------------

# Comprobar que hay genes significativos
sig_gene_names <- intersect(sig_genes$Gene, rownames(txi$counts))

if (length(sig_gene_names) >= 1) {
  
  counts_mat <- txi$counts[sig_gene_names, , drop = FALSE]
  
  # Crear annotation_col correctamente
  annotation_col <- data.frame(
    Condition = dds$Condition
  )
  rownames(annotation_col) <- colnames(counts_mat)
  
  pheatmap(
    counts_mat,
    cluster_rows = FALSE,   # pocos genes
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = annotation_col,
    main = "Heatmap of significant genes"
  )
  
} else {
  message("No significant genes available for heatmap.")
}

# ============================================================
# Boxplots genes significativos
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# Colores por condición 
cond_colors <- c(
  "Normopeso" = "lightblue",
  "Obeso2" = "coral"
)

# Crear boxplots para todos los genes
plots <- lapply(unique(counts_long$Gene), function(gene_name) {
  
  ggplot(counts_long %>% filter(Gene == gene_name),
         aes(x = Condition, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.colour = "black", outlier.size = 1) +
    stat_boxplot(geom = "errorbar", width = 0.25) +
    scale_fill_manual(values = cond_colors) +
    labs(
      title = paste("Expresión del gen", gene_name),
      x = "Condición",
      y = "Nivel de expresión génica"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(size = 7),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "gray50", fill = NA, linewidth = 1),
      plot.margin = unit(c(1,1,1,1), "lines")
    )
})

# Combinar todos los boxplots en una sola figura
box_overview <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title = "Expresión génica por condición",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )

box_overview

