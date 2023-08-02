# Combining and Contrasting Hardy-Weinberg Equilibrium  in Multi-Group Settings

This repository contains code related to our research aimed at extending the δ-based robust Hardy-Weinberg Equilibrium (HWE) test proposed by Zhang and Sun (2022) to combine and contrast single-nucleotide polymorphisms (SNPs) across groups.
Quality control in genome-wide association studies includes checking for Hardy-Weinberg Equilibrium (HWE), which poses a challenge in multi-group studies. Traditional methods, such as Pearson’s χ² test, offer limited interpretability without an effect size. This study proposes a robust allele-based HWE testing approach, leveraging data from 1604 individuals from the 1000 Genomes Project. The method facilitates effect size estimation and standard error for HWE tests. An inverse-variance meta-analysis and second-order omnibus meta-analysis were performed across four super-populations and two sexes. HWE comparison heterogeneity was assessed using Cochran’s Q. We identified 95 SNPs with genome-wide significant HWE departures (p-value<5e-8) considering both sex and population, with 26 of these showing significant heterogeneity. The results provide a more comprehensive perspective on HWE deviations across different populations and sexes.
## Data Analysis Framework
![data-pipeline](https://github.com/jaffaromain/1kgp-HWE-JR/assets/64853264/5f411423-3efc-4a00-a4b6-86372b57c847)

## Data
The data used in this study is obtained from the 1000 Genomes Project, which underwent prior quality control by The Centre of Applied Genomics. The data set comprises genotyping of individuals across 19 populations on the Illumina Omni2.5 chip. Further details regarding the data-cleaning process can be found in the scripts included in this repository.

### Data Required:
Data can be downloaded here: http://tcag.ca/tools/1000genomes.html

'sampleTable`: table of quality statistics per SNP 

'SnpTable': table of quality statistics per sample

'Indep': binary format files of independent samples

## Repository Files:
`reference_alleles.R`: Script for retrieving reference and alternative alleles. 
* Link to Ref. Allele data: https://www.well.ox.ac.uk/~wrayner/strand/RefAlt.html

`Functions-script.R` : Scripts for HWE and Meta-analyses functions.

`Population-HWE-Script`: Runs HWE separately for groups.

`Population-meta-analysis` : Runs Meta-analyses.

`results/` : Contains files with analysis results.


## Dependencies
PLINK v1.90

R 4.2.3 (or newer)

## Running the Analysis

*Clone the repository to your local machine

*Ensure you have all the required dependencies installed

*Navigate to the scripts/ directory and run scripts in the order reported above

## References
Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A., Bender, D., Maller,
J., Sklar, P., De Bakker, P. I., Daly, M. J., et al. (2007). Plink: a tool set for whole-
genome association and population-based linkage analyses. The American journal of
human genetics, 81(3):559–575.

R Core Team (2023). R: A Language and Environment for Statistical Computing. R
Foundation for Statistical Computing, Vienna, Austria.

Roslin, N. M., Weili, L., Paterson, A. D., and Strug, L. J. (2016). Quality control
analysis of the 1000 genomes project omni2. 5 genotypes. BioRxiv, page 078600.

Willer, C. J., Li, Y., and Abecasis, G. R. (2010). Metal: fast and efficient meta-analysis
of genomewide association scans. Bioinformatics, 26(17):2190–2191.

Zhang, L. and Sun, L. (2022a). A generalized robust allele-based genetic association
test. Biometrics, 78(2):487–498.

Zhang, L. and Sun, L. (2022b). Unifying genetic association tests via regression:
Prospective and retrospective, parametric and nonparametric, and genotype-and
allele-based tests. Canadian Journal of Statistics, 50(4):1321–1338.
