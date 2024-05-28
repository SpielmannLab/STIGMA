# STIGMA

## Introduction
NGS has made gene analysis routine in clinics, yet many variants in genes of unknown function are labeled as variants of uncertain significance (VUSs) and do not contribute to a diagnosis of rare diseases, until further validated by experimental verification. Gene prioritization can help narrow down the list of candidate genes under consideration. . 

Here we introduce single-cell tissue-specific gene prioritization using machine learning (STIGMA) to prioritize disease gene for congenital malformations. STIGMA predicts the disease-causing probability of genes based on their expression profiles across cell types, while considering the temporal dynamics during the embryogenesis of a healthy (wild-type) organism.
![alt text](https://github.com/SpielmannLab/STIGMA/blob/main/GraphicalAbstract.png?raw=true)


Following scripts were used collect the features, optimize and build the model. <br />

## Input Feature preprocessing<br />
### STEP1: Fetching Single Cell Features: <br />

1.	Load the Seurat object after clustering or subclustering. <br />
2.	Utilize the 'AverageExpression' function to acquire cluster-specific gene expression. <br />
3.	Employ the 'HVFInfo' function to compute variance in expression within each cluster. <br />
4.	Determine the percentage of cells expressing the gene in each sub-cluster (PrctCellExpringGene). <br />
5.	Calculate the fold-change in expression between each sub-cluster and the remaining cells (FoldChange). <br />
6.	Proceed with the trajectory analysis pipeline, following [Monocle's](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/) methodology for each partition/cluster. <br />
7.	Bin the generated pseudo time and incorporate it as metadata into the Seurat object. Then, compute the average expression of all genes for each pseudo time bin column in the metadata. This data serves as input to the featurePreprocess/Bsplines.R script (a sample dataset is provided for script input). <br />

### STEP2: Fetching Gene Intrinsic Properties: <br />

1.	Extract gene intrinsic properties, such as promoter GC content, by executing featurePreprocess/PromoterGC.R. <br />
2.	Obtain gene constraints, including metrics like pLI, pNull, pRec, syn_Z, mis_Z, and lof_Z, for protein-coding genes from gnomAD (v.2.1.1). <br />
3.	Retrieve gene GC content from Biomart. <br />

### STEP3: Creating the input matrix: <br />
1. Once the features are obtained, they are stored as a [tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/input.tsv) with column as features and rows as genes. <br />
2. Annotate the genes with positive and negative classes.


## STIGMA optimization and gene prioritization<br />
(1) STIGMA gene prediction model can be optimized by running model/RandomForest_optimization.py [input.tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/input.tsv)<br />
(2) STIGMA gene prediction model to predict test genes by running <br />
python model/rf_model.py --inputfeature_matrix=[input.tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/input.tsv) --candidategenes=[candidate_genes.tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/CandidateGene.tsv) --n_estimators=<Output from optimization> --max_depth=<Output from optimization> --min_samples_split=<Output from optimization> --min_samples_leaf=<Output from optimization> --max_features=<Output from optimization> --bootstrap=<Output from optimization> --n_neighbors=<Output from optimization> <br />

## STIGMA validation <br />
(1) To run explorative analysis based on Monarch Initiative validation/monarch_analysis.py <br />


[1] Absolute paths are used at certain instances, which will need to be adapted, as needed. <br />
[2] Anaconda was used to set up the necessary environments <br />

## Citations
1. Balachandran, S., Prada-Medina, C.A., Mensah, M.A., Kakar, N., Nagel, I., Pozojevic, J., Audain, E., Hitz, M.-P., Kircher, M., Sreenivasan, V.K.A., et al. (2024). STIGMA: Single-cell tissue-specific gene prioritization using machine learning. Am J Hum Genet, S0002-9297(23)00443-3. https://doi.org/10.1016/j.ajhg.2023.12.011.

