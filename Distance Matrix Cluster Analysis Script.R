# Load readxl
library(readxl)

# Read the distance matrix
dist_mat <- read_excel("Distance Matrix Phyml.xlsx")

# Convert the first column to row names
dist_mat <- as.data.frame(dist_mat)
rownames(dist_mat) <- dist_mat[[1]]
dist_mat[[1]] <- NULL

# Make sure it's numeric
dist_mat <- as.matrix(dist_mat)
mode(dist_mat) <- "numeric"

# Convert to a dist object
d <- as.dist(dist_mat)

# Perform hierarchical clustering (average linkage)
hc <- hclust(d, method = "average")

# Plot dendrogram
plot(hc, main = "Average Linkage Clustering", xlab = "", sub = "")

# Load metadata
metadata <- read_excel("metadata.xlsx")

# Make sure sample order matches distance matrix
metadata <- metadata[match(rownames(dist_mat), metadata$Sample), ]

# Extract true labels (e.g., "Site", "Cluster")
true_labels <- as.character(metadata$Cluster)

# Load mclust
library(mclust)

# Get cluster assignments from dendrogram
cluster_labels <- cutree(hc, k = length(unique(true_labels)))

# Calculate Adjusted Rand Index
rand_index <- adjustedRandIndex(cluster_labels, true_labels)

# Print result
cat("Adjusted Rand Index:", rand_index, "\n")

# RESULTS: Adjusted Rand Index: -0.08433735 

