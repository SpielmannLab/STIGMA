# STIGMA

Following scripts were used collect the features, optimize and build the model. <br />

## Input Feature preprocessing<br />
(1) The bsplines from the pseudotime bin data can be fit by running featurePreprocess/Bsplines.R <br />
(2) Gene interinsic property, promoter GC content can be mined by running featurePreprocess/PromoterGC.R <br />

## Model optimization and gene prioritization<br />
(1) STIGMA gene prediction model can be optimized by running model/RandomForest_optimization.py <br />
(2) STIGMA gene prediction model to predict test genes by running <br />
python model/rf_model.py --inputfeature_matrix=/Feature_input.tsv --candidategenes=candidate_genes.tsv --n_estimators=89 --max_depth=30 --min_samples_split=5 --min_samples_leaf=1 --max_features=auto --bootstrap=True --n_neighbors=5 <br />

## STIGMA validation <br />
(1) To run explorative analysis based on Monarch Initiative validation/monarch_analysis.py <br />


[1] Absolute paths are used at certain instances, which will need to be adapted, as needed. <br />
[2] Anaconda was used to set up the necessary environments <br />
