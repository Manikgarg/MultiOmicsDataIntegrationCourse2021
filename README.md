# Comparing different multi-omics data integration tools
 
## Overview
Several methods for multi-omics data integration exist, however, choosing the best method for a given dataset is still a challenge. In a recently available preprint, Laura et al.1 tried to benchmark nine different data integration methods and found that the methods Integrative NMF (intNMF)<sup>2<\sup> and Multiple co-inertia analysis (MCIA)3 performed consistently better than other seven methods across multiple datasets tested. However, these datasets had an equal number of samples in all modalities which is not usually the case. Moreover, the datasets tested did not account for count data or binary data and had no missing values. Therefore, the aim of the project will be to assess their performance on the chronic lymphocytic leukaemia (CLL)4 multi-omics dataset having unequal number of samples across different modalities, missing values and different omics data-types (such as binary mutation data as well as continuous RNA expression data) and compare its performance with the method, Multi-omics factor analysis (MOFA)5, specifically designed to account for these shortcomings.

## Dataset
CLL dataset contains information from 200 patients on
- somatic mutations
- RNA expression data
- DNA methylation
- ex vivo drug response

## Take-away
At the end of the course the participants will:
get hands-on experience in using atleast 3 different data-integration methods
learn advantages and limitations of different methods
get an intuition of which method to apply for which kind of dataset.
get hands-on experience on performing survival analysis, pathway analysis and clinical association analysis as a way to assess different methods.
get hands-on experience with handling missing values

## Requirements
Basic knowledge of R required.

## References
- Cantini, Laura, et al. “Benchmarking joint multi-omics dimensionality reduction approaches for cancer study.” bioRxiv(2020).
- Argelaguet, Ricard, et al. "Multi‐Omics Factor Analysis—a framework for unsupervised integration of multi‐omics data sets." Molecular systems biology 14.6 (2018).
- Chalise, Prabhakar, and Brooke L. Fridley. "Integrative clustering of multi-level ‘omic data based on non-negative matrix factorization algorithm." PloS one 12.5 (2017).
- Bady, Pierre, et al. "Multiple co-inertia analysis: a tool for assessing synchrony in the temporal variability of aquatic communities." Comptes rendus biologies 327.1 (2004): 29-36.
- Dietrich, Sascha, et al. "Drug-perturbation-based stratification of blood cancer." The Journal of clinical investigation 128.1 (2018): 427-445.
