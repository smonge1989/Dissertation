# ================================
# Required packages
# ================================
if (!require("readxl")) install.packages("readxl")
if (!require("mclust")) install.packages("mclust")
if (!require("pheatmap")) install.packages("pheatmap")

library(readxl)
library(mclust)
library(pheatmap)

# ================================
# Step 1: Load Distance Matrix
# ================================
dist_mat <- read_excel("Distance Matrix Phyml.xlsx")

dist_mat <- as.data.frame(dist_mat)
rownames(dist_mat) <- dist_mat[[1]]
dist_mat[[1]] <- NULL

dist_mat <- as.matrix(dist_mat)
mode(dist_mat) <- "numeric"

d <- as.dist(dist_mat)

# ================================
# Step 2: Perform Clustering
# ================================
hc <- hclust(d, method = "average")

plot(hc, main = "Average Linkage Clustering", xlab = "", sub = "")

# ================================
# Step 3: Load Metadata
# ================================
metadata <- read_excel("metadata.xlsx")

# Replace 'Sample' with actual column name that matches rownames
metadata <- metadata[match(rownames(dist_mat), metadata$Sample), ]

# Extract true labels (adjust column name if needed)
true_labels <- as.character(metadata$Cluster)

# ================================
# Step 4: Calculate Adjusted Rand Index
# ================================
k <- length(unique(true_labels))
cluster_labels <- cutree(hc, k = k)

if (any(is.na(true_labels))) {
  stop("NA values found in true labels. Check metadata.")
}

rand_index <- adjustedRandIndex(cluster_labels, true_labels)
cat("Adjusted Rand Index:", rand_index, "\n")

# ================================
# Step 5: Heatmap with Annotations
# ================================
# Build annotation dataframe for heatmap
annotation <- data.frame(Cluster = true_labels)
rownames(annotation) <- rownames(dist_mat)

annotation <- data.frame(Site = as.character(metadata$Cluster))
rownames(annotation) <- rownames(dist_mat)

# Use pheatmap
pheatmap(dist_mat,
         clustering_distance_rows = d,
         clustering_distance_cols = d,
         clustering_method = "average",
         annotation_row = annotation,
         annotation_col = annotation,
         main = "Heatmap with Cluster Annotations",
         border_color = NA)
