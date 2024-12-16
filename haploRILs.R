# Title: haploRILs.R
# Version: 1.0.0
# Author: jmontero
# Date: 2024-10-11
# Description: Founder haploblock reconstruction: (1) Converts heterozygous sites to missing calls, (2) Compares each progeny (RIL/DH) with all founders and creates score matrix with 0/1 scores per SNP and founder based on similarity to each founder, (3) Bin scores by windows with nSnp SNPs and step, (4) Detect highest-scoring founder, (5) Build haploblocks based on contiguous windows assigned to same founder, (6) Validate putative crossovers by comparing the highest-scoring founders between the CO with the K context regions before and after these, (7) Merge dissenting isolated haploblocks, (8) Return the haploblock pedigree information, score details and coordinates.
# Usage: ./haploRILs.R <filename> <nSnp> <step> <K> <number of founders>
# Input: PED/MAP called <filename>. The PED has individual by rows, founders first and then descendants, with ID in first col and SNP allele in the rest. haploRILs assumes the genotypes are represented in diploid format, with two numbers per loci, 1/2/0 for minor/major/missing allele (genotypes would be 1 1/1 2/2 2/0 0). Missing genotype = 0. The MAP file has SNPs by rows and chromosome and physical distance by columns. The number of cols in PED should be (number of rows in MAP file)*2 + 1
# Output: <filepath>_<nSnp>_<step>_<K>.hbk (haploblocks)
# Remarks: Now allows to set up the minimum number of SNPs per window based on the aim resolution. Step as input
# 24.10.28: Lines 120-121 - filename added as RIL

#### Load required packages, install if required ####
commonly_used_packages = c("data.table", "dplyr", "data.table", "ggplot2", "tidyr")
new_packages = c("stringr") # Introduce other packages to import
pkgs = c(commonly_used_packages, new_packages)
for (pkg in pkgs) { if (!suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE))) { install.packages(pkg, character.only = TRUE, repos = "https://cran.uni-muenster.de/", quiet = TRUE) ; require(pkg, character.only = TRUE, quietly = TRUE) } }

#### Parse argument from command line ####
args = commandArgs(trailingOnly = TRUE)
fileName = as.character(args[1])
nSnp = as.numeric(args[2])
step = as.numeric(args[3])
K = as.numeric(args[4])
numberFounders = as.numeric(args[5])
minimumDistanceMethod = as.character("null")

#### Load functions script ####
source("/vol/agcpgl/jmontero/RILs/Analyses/HaploblockReconstruction/haploRILs/scripts/haploRILs_function.R")

#### Define input and output names ####
inputName = paste0(fileName, ".ped")
outputName = paste0(fileName, "_", nSnp, "_", step, "_", K, ".hbk")

#### Import PED input and edit ####
ped = as.matrix(fread(inputName) %>% as.data.frame())
# Remove first col to leave only allele calls
rownames_ped = ped[, 1]
ped = ped[, -1]
# Make the data numeric
ped = apply(ped, 2, as.numeric)
# Add rownames
rownames(ped) = rownames_ped

#### Convert PED file format ####
# The input PED should be diploid and have the code
# 1/2/0 for minor/major/missing allele. Since we assume
# the samples must be mostly inbred, we can convert
# the genotyping format to haploid, with 1/2/0 denoting
# 1 1/2 2/0 0 in the previous format.
# Split the PED in two tables: first allele and second
# allele
ped1 = ped[,seq(1, ncol(ped), 2)] ; ped2 = ped[,seq(2, ncol(ped), 2)]
# Sum the numeric values of the tables
ped = ped1+ped2
# Convert the sum to the actually corresponding labels
ped = ifelse(ped == 2, 1, ifelse(ped == 4, 2, 0))

#### Import MAP file and edit ####
# Create SNP ID tags combining chromosome ID with SNP
# position (bp)
map = fread(paste0(fileName, ".map"), drop = c(2, 3)) %>%
  tidyr::unite("marker_id", c(V1, V4), sep = "_", remove = FALSE) %>%
  rename(chr = V1, pos = V4)
# Set variable chr
chr = map$chr[1]
# Add SNP ID as PED colnames
colnames(ped) = map$marker_id
# Save founders and RILs names into vectors
founderNames = rownames(ped)[1:numberFounders]
rilNames = rownames(ped)[(numberFounders+1):length(rownames(ped))]

#### Split the SNPs into windows ####
# Calculate density as Kilobasepairs per SNP
density_kbpsnp = ((map$pos[length(map$pos)]-map$pos[1])/1000)/nrow(map)
# Use density to calculate minimum distance per block
# proportional to the minimum number of SNPs (nSnps)
# depending on the minimumDistanceMethod = c(null,
# single, double, quadriple)
# (This is an old haploRILs functionality, should stay
# as "null")
minimum_distance_applied = ifelse(minimumDistanceMethod == "null", 0, ifelse(minimumDistanceMethod == "single", nSnp*density_kbpsnp*1000, ifelse(minimumDistanceMethod == "double", nSnp*density_kbpsnp*1000*2, ifelse(minimumDistanceMethod == "quadruple", nSnp*density_kbpsnp*1000*4))))
# The aim resolution is 2kbp. Determine the minimum_snps per windows that is necessary to reach that resolution 
#nSnp = defineMinimumNumberSNPs(map, aimResolution = 2e3, noLessThan_SNPs = 5)
# Define the coordinates of the windows where haplotypes
# will be assigned
slidingWindows = getSlidingWindows(map, minimum_snps = nSnp, minimum_distance = minimum_distance_applied, step = nSnp/step)

#### haploRILs core ####
haplorilsOutput = data.frame()
for (id in rilNames) {
  #### Assign scores to PED file ####
  idScores = assignScoresRules2(ped[c(founderNames, id),])
  #### Sum founder scores by sliding windows ####
  idSumScores = blockWindowBuilder(idScores,
                                   slidingWindows,
                                   corWindowLength = K,
                                   exportMatrix = TRUE)
  #### Extract a table with the most-likely founder per hblock, after filtering hblocks with the K system ####
  idBlocks = inferMostLikelyPathAlternative6(idSumScores,
                                             k = K,
                                             slidingWindows = slidingWindows,
                                             genotypes = idScores,
                                             export = "blocks")
  #### Finalize the table by adding context information ####
  idBlocks = idBlocks %>% 
    mutate(id = id,
    chr = chr,
    nSnp = nSnp,
    step = step,
    K = K) %>%
    select(id, chr, nSnp, step, K, everything())
  #### Filter blocks to merge repeating consecutive blocks and to eliminate isolated dissenting blocks ####
  filteredIdBlocks = filterBlocks(as.data.frame(idBlocks))
  haplorilsOutput = rbind(haplorilsOutput, filteredIdBlocks)
}

#### Add filename as RIL ####
haplorilsOutput = haplorilsOutput %>% dplyr::mutate(RIL = fileName) %>% select(RIL, everything())

#### Save output file ####
write.table(haplorilsOutput, outputName, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#### Print session info #####
sessionInfo()
