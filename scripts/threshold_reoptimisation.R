# Threshold re-optimisation for BCMA

library(ape)
library(adegenet)
library(spider)

## 1. Load your full, multi-specimen, delimited dataset
# Edit the file path below. Sequence names should follow <SpeciesCode>_<ID>,
# e.g. AA_MT895763, AB_NC068898, AL... so species labels can be extracted.
full_alignment <- fasta2DNAbin("data/full_delimited_dataset_COI.fasta", quiet = TRUE)

species_labels <- sapply(strsplit(dimnames(full_alignment)[[1]], "_"), `[`, 1)

## 2. Calculate pairwise distances
dist_matrix <- dist.dna(full_alignment, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
diag(dist_matrix) <- NA

## 3. Test a range of thresholds and find the one minimising cumulative error
threshVal <- seq(0.001, 0.02, by = 0.001)

sens <- lapply(threshVal, function(x) threshOpt(dist_matrix, species_labels, threshold = x))
sensMat <- do.call(rbind, sens)

# sensMat columns, in order: threshold, cumulative error, correct, incorrect, ambiguous
# (see ?spider::threshOpt for full column definitions)
barplot(t(sensMat)[4:5, ], names.arg = paste0(sensMat[, 1] * 100, "%"),
        main = "Incorrect (bottom) vs ambiguous (top) identifications per threshold")

new_threshold <- threshVal[which.min(sensMat[, 5])]
new_threshold

## 4. Check identification accuracy at the new threshold
bcma_check <- bestCloseMatch(dist_matrix, species_labels, threshold = new_threshold, names = TRUE)

# Report the percentage correctly assigned at this threshold, and inspect any
# incorrect or ambiguous cases individually before adopting the new value.
table(bcma_check[, "results"])  # column name may differ slightly by spider version - check first

## 5. Update the main script
# Once you are satisfied with the new threshold, manually update the relevant
# value (coi_threshold or nd2_threshold) in scripts/specimen_assignment.R,
# and update the "validated for" species list in this repository's README
# and in data/reference_accessions.csv to reflect the expanded reference set.
