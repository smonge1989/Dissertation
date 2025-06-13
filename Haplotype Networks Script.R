# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install necessary Bioconductor packages if not already installed
packages_needed <- c("ape", "pegas", "msa", "Biostrings", "GenomeInfoDbData", "GenomeInfoDb")
for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("msa", "Biostrings", "GenomeInfoDbData", "GenomeInfoDb")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Load libraries
library(ape)
library(pegas)
library(msa)
library(Biostrings)

# Set folder with FASTA files
folder_path <- "E:/DNAmt/DLoop"

# Create output folder for plots or aligned files
output_dir <- file.path(folder_path, "haplo_network_plots")
dir.create(output_dir, showWarnings = FALSE)

# Combine FASTA files into one multi-sequence FASTA
fasta_files <- list.files(folder_path, pattern = "\\.fasta$|\\.fa$", full.names = TRUE)
all_seqs <- DNAStringSet()

for (f in fasta_files) {
  seq <- readDNAStringSet(f)
  all_seqs <- append(all_seqs, seq)
}

# Save combined FASTA
combined_fasta <- file.path(folder_path, "dloop_combined_sequences.fasta")
writeXStringSet(all_seqs, filepath = combined_fasta)
cat("✅ Combined FASTA saved to:", combined_fasta, "\n")

# Align sequences using ClustalW
seqs <- readDNAStringSet(combined_fasta)
alignment <- msa(seqs, method = "ClustalW")

# Convert to DNAStringSet
alignment_seqinr <- msaConvert(alignment, type = "seqinr::alignment")
aligned_seqs <- alignment_seqinr$seq
aligned_names <- alignment_seqinr$nam
aligned <- DNAStringSet(aligned_seqs)
names(aligned) <- aligned_names

# Save aligned sequences
aligned_fasta <- file.path(folder_path, "dloop_aligned_sequences.fasta")
writeXStringSet(aligned, filepath = aligned_fasta)
cat("Alignment complete. Aligned file saved to:", aligned_fasta, "\n")

# Check alignment lengths
widths <- width(aligned)
print(widths)
if (length(unique(widths)) == 1) {
  cat("✅ All sequences are aligned and of equal length.\n")
} else {
  cat("❌ Sequences are NOT aligned properly (different lengths).\n")
}

# Clean sequence names for PopArt
names(aligned) <- gsub("[^A-Za-z0-9_]", "_", names(aligned))
cleaned_fasta <- file.path(folder_path, "dloop_aligned_sequences_cleaned.fasta")
writeXStringSet(aligned, filepath = cleaned_fasta)

# Convert to NEXUS
dna <- read.dna(aligned_fasta, format = "fasta")
rownames(dna) <- gsub("[^A-Za-z0-9_]", "_", rownames(dna))
output_nexus <- file.path(folder_path, "dloop_aligned_sequences.nex")
write.nexus.data(dna, file = output_nexus, format = "dna", interleaved = FALSE)
cat("✅ NEXUS file saved to:", output_nexus, "\n")

#Median-Joining Network constructed using PopArt, a separate program
