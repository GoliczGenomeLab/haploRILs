# haploRILs
> Version: 1.0.0
**Founder haplotype reconstruction with SNP data**
## Installation
```{shell Installation}
git clone https://github.com/GoliczGenomeLab/haploRILs.git
```
## Input
PED/MAP called `{filename}` installed in same dir as input.
- The PED file has individual by rows, founders first and then descendants, with ID in first col and SNP allele in the rest. haploRILs assumes the genotypes are represented in diploid format, with two numbers per loci, 1/2/0 for minor/major/missing allele (genotypes would be 1 1/1 2/2 2/0 0). Missing genotype = 0
- The MAP file has SNPs by rows and chromosome and physical distance by columns. The number of cols in PED should be `(number of rows in MAP file)*2 + 1`

## Usage
```{r Usage}
Rscript haploRILs.R {file name} {nSnp} {step} {K} {number of founders}
```
## Method
![combined_fill_1s](https://github.com/user-attachments/assets/916651a4-aed1-4702-8aca-a2631dcb83d6)
*The GIF shows two examples of founder haploblock reconstruction in chromosomes where a putative crossover is being detected between SNP windows 3 and 4 (K=1, step=1). On top, the crossover is confirmed because none of the pairs of SNP windows compared shares highest-scoring founders. Instead, the crossover below is rejected because the window after, 5 (i+k), shares founder with the windows before, 3 (i-1) and 2 (i-k-1). This context-based method for putative crossover validation aims to reduce the impact of genotyping errors or missing data on founder haplotype assignment.*

1. Heterozygous sites are converted to missing genotypes.
2. Descendants (RIL/DH) are compared with founders {PED lines 1 - `{number of founders}`} and a score matrix is created. The score matrix contains the values 0/1 scores per SNP and founder based on similarity to each founder.
3. Scores are summed by windows of `{nSnp}` SNPs and `{step}` overlap.
4. Highest-scoring founder are detected by windows.
5. Haploblocks are built based on contiguous windows assigned to same highest-scoring founder(s). Putative crossovers are marked in haploblock transitions.
	* Putative crossovers are validated by comparing the highest-scoring founders between the windows before and after the crossover with those of the `{K}`-context regions before and after these. If compared windows share highest-scoring founders, putative crossovers are rejected (watch GIFs).
6. Dissenting isolated haploblocks are filtered and imputed (AAABAAA --> AAAAAAA)
7. Return the haploblock pedigree information, score details and coordinates.

## Customizable parameters
|Parameter			|Description
|-------------------------------|---------------------------------------------------------------------------------------|
|nSnp                           |*Window size*, haploRILs splits the chromosome in windows of this number of SNPs	|
|step				|*Window overlap*, 1=Stable window, >1=Sliding window					|
|K                              |*Context window size*, number of SNPs on each side of the putative CO for validation	|
|number of founders		|Which line is the last founder in the PED file. The founders must go first in PED	|

## Output
One `.hbk` file containing founder haploblock information is generated in every haploRILs run named:
```{r Output}
{filepath}_{nSnp}_{step}_{K}.hbk
```
|Column name			|Explanation
|-------------------------------|-------------------------------------------------------------------------------|
|id				|Col 1 value in PED file							|
|chr				|Chromosome									|
|nSnp				|Window size used to build the hblock						|
|step				|Window step used to build the block. Equals the overlap of each window		|
|K				|Number of hblocks used for context-based validation of putative crossovers	|
|blocksFiltered			|Block index after filtering, so the last value is the number of retained blocks|
|block				|Block index before filtering, so it matches blocksFiltered if no filtering done|
|startIndex			|Hblock start position (line number in MAP file)				|
|endIndex			|Hblock end position (line number in MAP file)					|
|startBP			|Hblock start position (physical position of SNP in chr)			|
|endBP				|Hblock end position (physical position of SNP in chr)				|
|mostLikelyFoundersFiltered	|Highest-scoring founder after filtering					|
|mostLikelyFounders		|Highest-scoring founder before filtering					|
|mostLikelyFounders2		|2nd Highest-scoring founder before filtering					|
|blockSumScore1			|Block score of mostLikelyFounders						|
|blockSumScore2			|Block score of mostLikelyFounders2						|
|blockSumScoreRest		|Block score of mostLikelyFounders3						|
|startInterval			|Chromosome:{start bp position of first window}-Chromosome:{startBP, which is the end bp position of first window :raised_eyebrow:}
|endInterval			|Chromosome:{endBP, which is the start bp position of last window}-Chromosome:{end bp position of last window}|

An example of the output:
```{r Example output}
id      chr     nSnp    step    K       blocksFiltered  block   startIndex      endIndex        startBP endBP   mostLikelyFoundersFiltered      mostLikelyFounders  mostLikelyFounders2     blockSumScore1  blockSumScore2  blockSumScoreRest       startInterval   endInterval
123     chrSim  12      1       3       1       1       1       691     1       8281    Founder9    Founder9    Founder7    0.96562981186686    0.765197395079595       0.667020911722142       chrSim_1-chrSim_12      chrSim_8281-chrSim_8292
123     chrSim  12      1       3       2       2       692     1358    8304    16285   Founder7    Founder7    Founder4    0.863194002998501   0.809971814092954       0.696145512957807       chrSim_8293-chrSim_8304 chrSim_16285-chrSim_16296
123     chrSim  12      1       3       3       3       1359    1567    16308   18796   Founder8    Founder8    Founder14   0.839715311004785   0.832140669856459       0.611928537252222       chrSim_16297-chrSim_16308       chrSim_18793-chrSim_18796
```
# Notes

- `haploRILs.R` reads homemade functions from `haploRILs_function.R`. This is a well-documented file where functions can be tweaked to customize haploRILs default behaviour. It also contains unused functions.
- haploRILs is still in developmental stage so you are welcome to suggest changes and share your experience on the issues section.

# Contact
Author: Jose Antonio Montero Tena, PhD researcher

- [Email](jose.a.montero-tena@ab.uni-giessen.de)

- [GitHub](https://github.com/jamonterotena)
