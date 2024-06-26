# STIGMA

## Introduction
NGS has made gene analysis routine in clinics, yet many variants in genes of unknown function are labeled as variants of uncertain significance (VUSs) and do not contribute to a diagnosis of rare diseases, until further validated by experimental verification. Gene prioritization can help narrow down the list of candidate genes under consideration. . 

Here we introduce single-cell tissue-specific gene prioritization using machine learning (STIGMA) to prioritize disease gene for congenital malformations. STIGMA predicts the disease-causing probability of genes based on their expression profiles across cell types, while considering the temporal dynamics during the embryogenesis of a healthy (wild-type) organism.
![alt text](https://github.com/SpielmannLab/STIGMA/blob/main/GraphicalAbstract.png?raw=true)

## Dependencies
### R:
1. Seurat
2. Monocle3
3. dplyr
4. splines
5. biomaRt
6. AnnotationHub
7. GenomicRanges
8. BSgenome
9. tidyr
10. GenomicFeatures
11. TxDb.Hsapiens.UCSC.hg38.knownGene

### Python
1. sklearn
2. imblern
3. pandas
4. numpy
5. math
6. sys
7. itertools
8. matplotlib

Inorder to run STIGMA we need to prepare the input matrix which consists of single cell features and gene intrinsic propersties, which can be done by following **STEP1** and **STEP2**. **STEP3** Formats the input. <br />


## Input Feature preprocessing<br />
### STEP1a: Fetching Single Cell Features from Seurat: <br />

1.	Load the [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) object after clustering or subclustering. <br />
2.	Utilize the ['AverageExpression'](https://satijalab.org/seurat/reference/averageexpression) function to acquire cluster-specific gene expression. <br />
3.	Employ the ['HVFInfo'](https://satijalab.org/seurat/reference/hvfinfo.sctassay) function to compute variance in expression within each cluster. <br />
4.	Determine the percentage of cells expressing the gene in each sub-cluster ([PrctCellExpringGene](https://rdrr.io/github/vertesy/Seurat.utils/man/PrctCellExpringGene.html)). <br />
5.	Calculate the fold-change in expression between each sub-cluster and the remaining cells ([FoldChange](https://satijalab.org/seurat/reference/foldchange)). <br />


### STEP1b: Fetching Single Cell temporal Features from Monocle3: <br />
1.	Proceed with the trajectory analysis pipeline, following [Monocle's](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/) methodology for each partition/cluster. <br />
2.	Bin the generated pseudo time and incorporate it as metadata into the Seurat object.
   ````

  pseudo <- data.frame(pseudotime(monocleobject, reduction_method = "UMAP"))
  colnames(pseudo)[1] <- 'pseudotime'
  seuratobject@meta.data$pseudotime<-pseudo$pseudotime[match(colnames(monocleobject), rownames(seuratobject))]
  seuratobject@meta.data$pseudotime.bin <- as.numeric(findInterval(seuratobject$pseudotime, quantile(seuratobject$pseudotime, seq(0,1, 1/10),na.rm=T)))

  ````
4.	Then, compute the average expression of all genes for each pseudo time bin column in the metadata.
````

avg_expr <- data.frame(AverageExpression(object=seuratobject, assays='RNA', slot='counts', group.by='pseudotime.bin'))
write.table(avg_expr,'Input_bsplines.tsv', sep='\t')

````
5.	Fit a spline by parsing [Input_bsplines.tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/Input_bsplines.tsv) <br />
   ````

   Rscript featurePreprocess/Bsplines.R Input_bsplines.tsv. 

````

### STEP2: Fetching Gene Intrinsic Properties: <br />

1.	Extract gene intrinsic properties, such as promoter GC content, by executing <br />
````
  	Rscript featurePreprocess/PromoterGC.R. 
````
3.	Obtain gene constraints, including metrics like pLI, pNull, pRec, syn_Z, mis_Z, and lof_Z, for protein-coding genes from [gnomAD](https://gnomad.broadinstitute.org/downloads#v4-constraint). In the paper we have used an older version of gnomad(v.2.1.1). <br />
4.	Retrieve gene GC content from [Biomart](https://www.ensembl.org/biomart/martview). <br />

### STEP3: Creating the input matrix: <br />
1. Once the features are obtained, they are stored as a [tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/input.tsv) with column as features and rows as genes. <br />
2. Annotate the genes with positive([Known disease gene](https://panelapp.genomicsengland.co.uk/panels/384/)) and negative classes(House keeping genes).


## STEP4: STIGMA optimization and gene prioritization<br />
1. STIGMA gene prediction model can be optimized by parsing [input.tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/input.tsv) <br />
   ````
   python3 model/RandomForest_optimization.py input.tsv
````
2. STIGMA gene prediction model to predict test genes with [input.tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/input.tsv) by running <br />
````
python3 model/rf_model.py --inputfeature_matrix=input.tsv --candidategenes=[candidate_genes.tsv](https://github.com/SpielmannLab/STIGMA/blob/main/sample_dataset/CandidateGene.tsv) --n_estimators=\<Output from optimization\> --max_depth=\<Output from optimization\> --min_samples_split=\<Output from optimization\> --min_samples_leaf=\<Output from optimization\> --max_features=\<Output from optimization\> --bootstrap=\<Output from optimization\> --n_neighbors=\<Output from optimization\> <br />
````
## STEP5: STIGMA validation <br />
1. To optain the phenotypes associated with the congenital disease based on Monarch Initiative <br />
````
python3 validation/monarch_analysis.py 
````

[1] Absolute paths are used at certain instances, which will need to be adapted, as needed. <br />
[2] Anaconda was used to set up the necessary environments <br />

## Citations
1. Balachandran, S., Prada-Medina, C.A., Mensah, M.A., Kakar, N., Nagel, I., Pozojevic, J., Audain, E., Hitz, M.-P., Kircher, M., Sreenivasan, V.K.A., et al. (2024). STIGMA: Single-cell tissue-specific gene prioritization using machine learning. Am J Hum Genet, S0002-9297(23)00443-3. https://doi.org/10.1016/j.ajhg.2023.12.011.

