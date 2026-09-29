# Acroteriobatus specimen assignment tool

This tool accompanies:

> Groeneveld, M.J., Klein, J.D., van Staden, M., Bennett, R.H., Dicken, M.L.,
> Ebert, D.A., Mann, B.Q., Perumal, K., Watson, R.G.A. & Bester-van der Merwe, A.E.
> (2026). Molecular taxonomy of the guitarfish genus *Acroteriobatus*
> (Rhinopristiformes: Rhinobatidae). *Marine Biodiversity*.
> https://doi.org/10.1007/s12526-025-01617-x

If you use this tool, please cite the paper above. An example of the
consensus-assignment output is given in Tables S5 and S6 of the paper's
supplementary information.

## What this does

Given fasta alignments of unidentified specimens (COI and/or ND2), this
workflow assigns each one to a described *Acroteriobatus* species using four
methods:

| Method | Package | Approach |
|---|---|---|
| BCMA | `spider` | Best Close Match, at a fixed, pre-validated distance threshold |
| BP | `BarcodingR` | Probabilistic (Bayesian-like) assignment |
| FZ | `BarcodingR` | Fuzzy-set identification |
| FZKMER | `BarcodingR` | Fuzzy-set identification using k-mers |

Each method is run separately and its results are written to its own file.
Groeneveld et al. (2026) treat a specimen as robustly assigned where at
least 3 of the 4 methods agree, and as ambiguous where they do not - see the
paper for how these were compared and combined. Ambiguous cases should in any case be
checked against independent evidence (e.g. morphology) rather than forced
to a single species.

Species delimitation (i.e. defining the MOTUs/species groups themselves) is
**not** performed by this script. This tool assumes species boundaries are
already established and assigns new, unidentified sequences against
that reference framework.

## Important: BCMA uses the full alignment; BP, FZ and FZKMER use haplotypes

Your unidentified specimens need to be prepared as **two separate files per
marker**:

- **A full alignment**, containing every individual specimen, used for BCMA.
  BCMA works directly on pairwise genetic distances between individual
  sequences, so it is run on the complete dataset.
- **A haplotype alignment**, containing one representative sequence per
  unique haplotype to avoid redundant computation, used for BP, FZ and FZKMER.
  Once you have the result for a given haplotype, it applies to every specimen that
  shares it

## Important: how the BCMA threshold works here

The BCMA distance thresholds used by `specimen_assignment.R` (0.4% for COI,
1.1% for ND2) are fixed values, taken directly from the validation carried
out in Groeneveld et al. (2026). They are **not** re-calculated each time
you run the script.

The reference set within this tool contains only one sequence per species,
and threshold optimisation needs several known specimens per species to
test against (it works by checking whether each specimen's nearest
neighbour is of the same species - with only one sequence per species, that
comparison doesn't exist). Threshold optimisation therefore has to be done
separately, on a full, multi-specimen dataset, not on the one-per-species
reference set used for routine assignment - see
`scripts/threshold_reoptimisation.R`.

**The fixed thresholds are validated for a specific species set, not all
eight described species:**

- COI: validated for *A. annulatus*, *A. blochii*, *A. andysabini*,
  *A. leucospilus*, *A. variegatus* and *A. zanzibarensis*. No COI sequences
  currently exist for *A. omanensis* or *A. salalah*, so the COI threshold
  has not been tested against these two species.
- ND2: validated across all eight currently sequenced species (i.e. all of
  the above plus *A. omanensis* and *A. salalah*).
- Neither marker's reference set currently includes *A. ocellatus* or
  *A. stehmanni* (no sequence data available for either species at time of
  writing).

If sequences for any of these missing species become available and are
added to the reference set, re-run `scripts/threshold_reoptimisation.R`
before trusting the assignment results for the expanded set - I will also
update this script as new sequences become available.

## Repository structure

This repository contains only the scripts and the reference data:

```
.
├── README.md
├── scripts/
│   ├── specimen_assignment.R          # main workflow - run this for routine assignment
│   └── threshold_reoptimisation.R     # only run if the reference set changes
└── data/
    ├── reference_sequences_COI.fasta
    ├── reference_sequences_ND2.fasta
    └── reference_accessions.csv           # species <-> GenBank accession lookup
```

Create and store your own unidentified specimens and the results the script produces
locally in this structure:

```
data/
├── unknown_COI_full.fasta          # every specimen - for BCMA
├── unknown_COI_haplotypes.fasta    # one sequence per haplotype - for BP, FZ, FZKMER
├── unknown_ND2_full.fasta
└── unknown_ND2_haplotypes.fasta

output/                             # created automatically when you run the script
```

## Getting the reference sequences

Two ways to obtain the reference fasta files:

1. **Download directly from this repository** - `reference_sequences_COI.fasta`
   and `reference_sequences_ND2.fasta` in `data/` contain one representative
   sequence per described species used in Groeneveld et al. (2026).
2. **Download from GenBank/BOLD yourself**, using the accession numbers
   listed in `data/reference_accessions.csv` (all other accessions given in Table S2
   of the paper's supplementary information). In-house generated sequences are
   additionally available under project **MTACR** on BOLD:
   https://portal.boldsystems.org/result?query=MTACR[recordsetcode]

Sequence names in the reference fasta files follow the convention
`<SpeciesCode>_<AccessionNumber>` (e.g. `AA_PV814390` for *A. annulatus* COI),
matching the `Species code` column in `reference_accessions.csv`.

## Requirements

R (>= 4.0) and the following packages:

```r
install.packages(c("ape", "adegenet", "dplyr"))
install.packages("spider", repos = "http://R-Forge.R-project.org")
```

`BarcodingR` is not on CRAN - install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("qinguoyou/BarcodingR")
```

## Usage

1. Clone this repository structure.
2. Place the reference fasta files in `data/`.
3. Prepare your own unidentified specimens locally as fasta alignments
   and place them in `data/`.
4. Run `scripts/specimen_assignment.R`. It processes COI and ND2 as two
   separate sections in the same script - if you only have
   data for one marker, just run that section.
5. Each method's results are written to its own file in `output/`,
   created automatically when you run the script.

## Contact

Questions or issues: open a GitHub issue on this repository, or contact me
on ResearchGate: https://www.researchgate.net/profile/Mia-Groeneveld.
