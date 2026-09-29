## Acroteriobatus specimen assignment tool
# IMPORTANT - read the README before running

## Required packages
# install.packages(c("ape", "adegenet", "dplyr"))
# install.packages("spider", repos = "http://R-Forge.R-project.org")
# BarcodingR is not on CRAN - install from GitHub:
# devtools::install_github("qinguoyou/BarcodingR")

library(ape)         # dist.dna()
library(adegenet)    # fasta2DNAbin()
library(spider)       # bestCloseMatch()
library(BarcodingR)   # barcoding.spe.identify(), barcoding.spe.identify2()

## PART A: COI

## --- A1. File paths (edit if your file names/locations differ) ------------
coi_reference_fasta <- "data/reference_sequences_COI.fasta"
coi_query_full_fasta <- "data/unknown_COI_full.fasta"          # every specimen - for BCMA
coi_query_hap_fasta  <- "data/unknown_COI_haplotypes.fasta"    # one per haplotype - for BP/FZ/FZKMER

## --- A2. Fixed, validated threshold -----------------------------------------
coi_threshold <- 0.004   # 0.4%, from Groeneveld et al. (2026)

## --- A3. Load sequences ------------------------------------------------------
coi_reference  <- fasta2DNAbin(coi_reference_fasta, quiet = TRUE)
coi_query_full <- fasta2DNAbin(coi_query_full_fasta, quiet = TRUE)
coi_query_hap  <- fasta2DNAbin(coi_query_hap_fasta, quiet = TRUE)

# Sequence names should follow the convention <SpeciesCode>_<AccessionNumber>
# for reference sequences (e.g. AA_MT895763), and any consistent naming for
# query sequences (e.g. unknown_001).

## --- A4. Method 1: BCMA (spider), run on the FULL alignment ----------------
coi_combined_full <- rbind(coi_reference, coi_query_full)

coi_ref_labels        <- sapply(strsplit(dimnames(coi_reference)[[1]], "_"), `[`, 1)
coi_query_full_labels <- rep("unknown", nrow(coi_query_full))
coi_full_labels        <- c(coi_ref_labels, coi_query_full_labels)

coi_dist <- dist.dna(coi_combined_full, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
diag(coi_dist) <- NA

coi_bcma_all <- bestCloseMatch(coi_dist, coi_full_labels, threshold = coi_threshold, names = TRUE)

# Keep only the query rows - reference rows were already validated in the paper.
coi_bcma <- coi_bcma_all[coi_full_labels == "unknown", ]

write.csv(coi_bcma, "output/coi_bcma_results.csv", row.names = FALSE)

## --- A5. Methods 2-4: BP, FZ, FZKMER (BarcodingR), run on the HAPLOTYPE alignment ---
coi_bp <- barcoding.spe.identify(coi_reference, coi_query_hap, method = "bpNewTraining")
write.csv(coi_bp, "output/coi_bp_results.csv", row.names = FALSE)

coi_fz <- barcoding.spe.identify(coi_reference, coi_query_hap, method = "fuzzyId")
write.csv(coi_fz, "output/coi_fz_results.csv", row.names = FALSE)

coi_fzkmer <- barcoding.spe.identify2(coi_reference, coi_query_hap, kmer = 5, optimization = TRUE)
write.csv(coi_fzkmer, "output/coi_fzkmer_results.csv", row.names = FALSE)

## PART B: ND2

## --- B1. File paths ----------------------------------------------------------
nd2_reference_fasta <- "data/reference_sequences_ND2.fasta"
nd2_query_full_fasta <- "data/unknown_ND2_full.fasta"
nd2_query_hap_fasta  <- "data/unknown_ND2_haplotypes.fasta"

## --- B2. Fixed, validated threshold -----------------------------------------
nd2_threshold <- 0.011   # 1.1%, from Groeneveld et al. (2026)

## --- B3. Load sequences -------------------------------------------------------
nd2_reference  <- fasta2DNAbin(nd2_reference_fasta, quiet = TRUE)
nd2_query_full <- fasta2DNAbin(nd2_query_full_fasta, quiet = TRUE)
nd2_query_hap  <- fasta2DNAbin(nd2_query_hap_fasta, quiet = TRUE)

## --- B4. Method 1: BCMA (spider), run on the FULL alignment ----------------
nd2_combined_full <- rbind(nd2_reference, nd2_query_full)

nd2_ref_labels        <- sapply(strsplit(dimnames(nd2_reference)[[1]], "_"), `[`, 1)
nd2_query_full_labels <- rep("unknown", nrow(nd2_query_full))
nd2_full_labels        <- c(nd2_ref_labels, nd2_query_full_labels)

nd2_dist <- dist.dna(nd2_combined_full, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
diag(nd2_dist) <- NA

nd2_bcma_all <- bestCloseMatch(nd2_dist, nd2_full_labels, threshold = nd2_threshold, names = TRUE)
nd2_bcma <- nd2_bcma_all[nd2_full_labels == "unknown", ]

write.csv(nd2_bcma, "output/nd2_bcma_results.csv", row.names = FALSE)

## --- B5. Methods 2-4: BP, FZ, FZKMER (BarcodingR), run on the HAPLOTYPE alignment ---
nd2_bp <- barcoding.spe.identify(nd2_reference, nd2_query_hap, method = "bpNewTraining")
write.csv(nd2_bp, "output/nd2_bp_results.csv", row.names = FALSE)

nd2_fz <- barcoding.spe.identify(nd2_reference, nd2_query_hap, method = "fuzzyId")
write.csv(nd2_fz, "output/nd2_fz_results.csv", row.names = FALSE)

nd2_fzkmer <- barcoding.spe.identify2(nd2_reference, nd2_query_hap, kmer = 5, optimization = TRUE)
write.csv(nd2_fzkmer, "output/nd2_fzkmer_results.csv", row.names = FALSE)
