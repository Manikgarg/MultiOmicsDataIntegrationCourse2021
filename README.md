# Comparing different multi-omics data integration tools

## Overview
Several methods for multi-omics data integration exist, however, choosing the best method for a given dataset is still a challenge. In a recently available preprint, Laura et al.<sup>1</sup> tried to benchmark nine different data integration methods on TCGA cancer data and found that the methods Multi-omics factor analysis (MOFA)<sup>2</sup>, Multiple co-inertia analysis (MCIA)<sup>3</sup>, Joint and Individual Variation Explained (JIVE)<sup>4</sup> and Regularized Generalized Canonical Correlation Analysis (RGCCA)<sup>5</sup> performed consistently better than other five methods. However, this benchmark required matching of samples across different modalities, did not account for their performance on binary count data and did not account for the variations introduced by different imputation methods. Therefore, the aim of this project will be to benchmark the first three multi-omics data-integration methods on the chronic lymphocytic leukaemia (CLL)<sup>6</sup> dataset.

## Dataset
CLL dataset contains information from 200 patients on
- somatic mutations
- RNA expression data
- DNA methylation
- ex vivo drug response

## Take-away
At the end of the course the participants will:
- get hands-on experience in using three different data-integration methods
- learn advantages and limitations of different methods
- get an intuition of which method to apply for which kind of dataset.
- get hands-on experience with handling missing values

## Requirements
Working knowledge of R, basic understanding of maths/statistics and familiarity with gene-set enrichment analysis and survival analysis required. 

## Pre-course homework
Please read reference 1 to better understand the course. Going through references 2, 3 and 4 would be great but not required.

## References
1. Cantini, Laura, et al. “Benchmarking joint multi-omics dimensionality reduction approaches for cancer study.” bioRxiv(2020).
2. Argelaguet, Ricard, et al. "Multi‐Omics Factor Analysis—a framework for unsupervised integration of multi‐omics data sets." Molecular systems biology 14.6 (2018).
3. Bady, Pierre, et al. "Multiple co-inertia analysis: a tool for assessing synchrony in the temporal variability of aquatic communities." Comptes rendus biologies 327.1 (2004): 29-36.
4. Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.
5. Tenenhaus, Arthur, and Michel Tenenhaus. "Regularized generalized canonical correlation analysis." Psychometrika76.2 (2011): 257.
6. Dietrich, Sascha, et al. "Drug-perturbation-based stratification of blood cancer." The Journal of clinical investigation 128.1 (2018): 427-445.
