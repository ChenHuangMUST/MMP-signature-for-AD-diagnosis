# MMP-signature-for-AD-diagnosis

## Description
Alzheimer's disease (AD) is an incurable neurodegenerative disorder. Early screening, particularly in blood has been regarded as an effective approach for AD diagnosis and prevention. In addition, metabolic dysfunction has been demonstrated to be closely related to AD, which might be reflected in whole blood transcriptome. Hence, we hypothesized that establishment of diagnostic model based on metabolic signatures of blood is a workable strategy. To the end, we initially constructed metabolic pathway pairwise (MPP) signatures to characterize the interplay among metabolic pathways. Then, a series of bioinformatic methodologies, e.g., differential expression analysis, functional enrichment analysis, and network analysis, etc. were used to investigate the molecular mechanism behind AD. Moreover, an unsupervised clustering analysis based on the MPP signature profile via Non-Negative Matrix Factorization (NMF) algorithm was utilized to stratify AD pa-tients. Finally, aimed at distinguishing AD patients from non-AD group, a metabolic path-way-pairwise scoring system (MPPSS) was established using multi-machine learning methods. As a result, many metabolic pathways correlated to AD were disclosed, including oxidative phos-phorylation and fatty acid biosynthesis, etc. NMF clustering analysis divided AD patients into two subgroups (S1 and S2), which exhibit distinct activities of metabolism and immunity. Typically, oxidative phosphorylation in S2 exhibits a lower activity than that in S1 and non-AD group, suggesting  patients in S2 might possess a more compromised brain metabolism. Additionally, Immune infiltration analysis showed that the patients in S2 might have phenomena of immune suppression, compared with S1 and non-AD group. These findings indicated that S2 probably have a more severe progression of AD. Finally, MPPSS could achieve an AUC of 0.73 in training dataset, 0.71 in testing dataset and an AUC of 0.82 on weighted average in five external validation datasets. Overall, our study successfully established a novel metabolism-based scoring system for AD di-agnosis using blood transcriptome, and provided new insight into the molecular mechanism of metabolic dysfunction implicated in AD.

## MPPSS
Source codes for generating results of "Metabolic Pathway Pairwise-based signature as a potential non-invasive diagnostic marker in AD patients"



## Requirements
- Boruta (8.0.0)
- xgboost (1.6.0.1)
- CIBERSORT (0.1.0)
- GSVA (1.44.5)
- ggfortify (0.4.15)
- ComplexHeatmap (2.12.1)
- preprocessCore (1.58.0)
- e1071 (1.7-12)
- randomForest (4.7-1.1)
- org.Hs.eg.db (3.15.0)
- AnnotationDbi (1.60.0)
- IRanges (2.32.0)
- S4Vectors (0.36.0)
- Biobase (2.58.0)
- BiocGenerics (0.44.0)
- clusterProfiler (4.4.4)
- pheatmap (1.0.12)
- pROC (1.18.0)
- ggpubr (0.4.0)
- ROCR (1.0-11)
- survival (3.3-1)
- glmnet (4.1-4)
- Matrix (1.5-3)
- caret (6.0-93)
- lattice (0.20-45)
- ggplot2 (3.4.0)




## Installation
All packages can be installed via bioconductor (https://www.bioconductor.org/) and CRAN. Generally, a couple of minutes is needed for installing each package.


## MPPSS predictions
- Code for reproducing multi machine learning methods and predictions are provided under the 'ML' folder.
- The expected run time is needed for running multi machine learning.
- NMF clustering and enrichment annotation are provided under the main folder.
