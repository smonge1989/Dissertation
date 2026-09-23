# Load required packages
library(readxl)
library(ape)        # for pcoa
library(ggplot2)

# Step 1: Read the distance matrix
dist_mat <- read_excel("Distance Matrix Phyml.xlsx")

# Step 2: Clean and format
dist_mat <- as.data.frame(dist_mat)
rownames(dist_mat) <- dist_mat[[1]]
dist_mat[[1]] <- NULL
dist_mat <- as.matrix(dist_mat)
mode(dist_mat) <- "numeric"

# Step 3: Convert to dist object
d <- as.dist(dist_mat)

# Step 4: Perform PCoA
pcoa_result <- pcoa(d)

# Step 5: Extract scores
scores <- as.data.frame(pcoa_result$vectors)
scores$Sample <- rownames(scores)

# Step 6: Plot the first two axes
ggplot(scores, aes(x = Axis.1, y = Axis.2, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5) +
  xlab(paste0("PCoA 1 (", round(pcoa_result$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste0("PCoA 2 (", round(pcoa_result$values$Relative_eig[2]*100, 1), "%)")) +
  ggtitle("PCoA based on Distance Matrix") +
  theme_minimal()

# Load metadata
metadata <- read_excel("metadata.xlsx")

# Merge metadata with PCoA scores
scores <- merge(scores, metadata, by.x = "Sample", by.y = "Sample")

# Plot with color by Site
ggplot(scores, aes(x = Axis.1, y = Axis.2, color = Site, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5) +
  xlab(paste0("PCoA 1 (", round(pcoa_result$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste0("PCoA 2 (", round(pcoa_result$values$Relative_eig[2]*100, 1), "%)")) +
  ggtitle("PCoA Colored by Site") +
  theme_minimal()



# Plot PCoA with ellipses around sites
ggplot(scores, aes(x = Axis.1, y = Axis.2, color = Site, fill = Site)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, color = NA) +  # Filled ellipses
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -0.5, hjust = 0.5) +
  xlab(paste0("PCoA 1 (", round(pcoa_result$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste0("PCoA 2 (", round(pcoa_result$values$Relative_eig[2]*100, 1), "%)")) +
  ggtitle("PCoA with Group Ellipses by Site") +
  theme_minimal()
