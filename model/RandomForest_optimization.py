import sys
import math
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_score,auc, recall_score, roc_auc_score, roc_curve
import itertools
from sklearn.preprocessing import MinMaxScaler
from imblearn.ensemble import BalancedBaggingClassifier
from sklearn.model_selection import GridSearchCV
from imblearn.over_sampling import BorderlineSMOTE, KMeansSMOTE, ADASYN, SMOTE
from sklearn.ensemble import VotingClassifier
from sklearn.impute import IterativeImputer
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from imblearn.pipeline import Pipeline, make_pipeline

def main():
	feat=sys.argv[1]

	df = pd.read_csv(feat,sep='\t')

	df = df[df['hg_ensembl_id'].notna()]
	df = df.set_index(['hg_ensembl_id','gene_symbol'])

	df1= df[df['Disease_Association'].isna()]

	X_test=df1.drop(['Disease_Association'],axis=1)


	df = df[df['Disease_Association'].notna()]
	plt.figure(figsize=(8, 8))
	sns.countplot('Disease_Association', data=df)
	plt.title('Classes Distribution')

	plt.savefig('ClassesDistribution.pdf',bbox_inches='tight')
	plt.clf()

	factor = pd.factorize(df['Disease_Association'], sort=True)
	df.Disease_Association = factor[0]
	definitions = factor[1]
	print(df.Disease_Association.head())
	print(definitions)
	Y = df[['Disease_Association']]
	X = df.drop(['Disease_Association'],axis=1)





	X.to_csv('train.tsv', sep='\t')
    # Impute missing values
	imp=IterativeImputer(max_iter=10, random_state=0)
	X_values=imp.fit_transform(X.values)
	X_test_values=imp.transform(X_test.values)
    #Scale the data
	scaler= MinMaxScaler()
	X_values=scaler.fit_transform(X_values)
	X_test_values=scaler.transform(X_test_values)

	X=pd.DataFrame(X_values,index=X.index, columns=X.columns)
	X_test=pd.DataFrame(X_test_values, index=X_test.index, columns=X_test.columns)

	# Number of trees in random forest
	n_estimators = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200],
	# Number of features to consider at every split
	max_features = ['auto', 'sqrt']
	# Maximum number of levels in tree
	max_depth = [10, 15, 20, 30, 40, 50, 75, 100, 150, 200, 250, 300, None]
	max_depth.append(None)
	# Minimum number of samples required to split a node
	min_samples_split = [2, 5, 10]
	# Minimum number of samples required at each leaf node
	min_samples_leaf = [1, 2, 4]
	# Method of selecting samples for training each tree
	bootstrap = [True, False]
    # Number of neigbors to consider to simulate a new datapoint
    n_neighbor = [2, 5, 10]

	new_params = {'randomforestclassifier__n_estimators': n_estimators,
               'randomforestclassifier__max_features': max_features,
               'randomforestclassifier__max_depth': max_depth,
               'randomforestclassifier__min_samples_split': min_samples_split,
               'randomforestclassifier__min_samples_leaf': min_samples_leaf,
               'randomforestclassifier__bootstrap': bootstrap,
               'adasyn__n_neighbors': n_neighbor}


	imba_pipeline = make_pipeline(ADASYN(random_state=42),RandomForestClassifier())


	rf_random = GridSearchCV(imba_pipeline, param_grid = new_params, cv = 3,scoring='recall')

	search=rf_random.fit(X, Y.values.ravel())
	print(search.best_params_)

if __name__ == "__main__":
	main()
