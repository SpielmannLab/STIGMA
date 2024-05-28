import argparse, sys
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LassoCV
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import precision_score,auc, recall_score, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import itertools
import os
from numpy import argmax
from sklearn.preprocessing import MinMaxScaler
from imblearn.over_sampling import BorderlineSMOTE, KMeansSMOTE, ADASYN, SMOTE
from sklearn.ensemble import VotingClassifier
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
def main():
	parser=argparse.ArgumentParser()
	parser.add_argument("--inputfeature_matrix", help="inputfeature_matrix with all the features from single cell and database as a .tsv file")
	parser.add_argument("--candidategenes", help="candidategenes as a tsv file")
	parser.add_argument("--n_estimators", help="n_estimator output from the optimization script")
	parser.add_argument("--max_depth", help="max_depth output from the optimization script")
	parser.add_argument("--min_samples_split", help="min_samples_split output from the optimization script")
	parser.add_argument("--min_samples_leaf", help="min_samples_leaf output from the optimization script")
	parser.add_argument("--max_features", help="max_features output from the optimization script")
	parser.add_argument("--bootstrap", help="bootstrap output from the optimization script")
	parser.add_argument("--n_neighbors", help="n_neighbors output from the optimization script")
	global args
	args=parser.parse_args()
	feat=args.inputfeature_matrix
	candidate_genes=args.candidategenes
	candidate_genes_df = pd.read_csv(candidate_genes,sep='\t')
	df = pd.read_csv(feat,sep='\t')
	print(df.columns)
	df = df[df['ensembl_id'].notna()]
	df = df.set_index(['ensembl_id','gene_symbol'])


	X_test=df.drop(['Disease_Association'],axis=1)


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
	print(f"Args: {args}\nCommand Line: {sys.argv}\nfoo: {args.n_neighbors}")
	Y = df[['Disease_Association']]
	X = df.drop(['Disease_Association'],axis=1)

	X = X.fillna(0)
	X_test = X_test.fillna(0)


	X.to_csv('train.tsv', sep='\t')

	imp=IterativeImputer(max_iter=10, random_state=0)
	X_values=imp.fit_transform(X.values)
	X_test_values=imp.transform(X_test.values)

	scaler= MinMaxScaler()
	X_values=scaler.fit_transform(X_values)
	X_test_values=scaler.transform(X_test)

	X=pd.DataFrame(X_values,index=X.index, columns=X.columns)
	X_test=pd.DataFrame(X_test_values, index=X_test.index, columns=X_test.columns)
	plt.figure(figsize=(12,10))
	cor = df.corr()
	cor.to_csv('Feature_correlation.tsv',sep='\t')
	cor_target = abs(cor["Disease_Association"])
	relevant_features = cor_target[cor_target>0.5]
	rel_cor = cor.loc[relevant_features.index, relevant_features.index]
	sns.heatmap(rel_cor, annot=True, cmap=plt.cm.Reds)


	kf = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)
	counter=0
	prediction=[]
	prediction_prob_0=[]
	prediction_prob_1=[]

	validation=[]
	validation_prob_0=[]
	validation_prob_1=[]

	validate_gene=[]
	validate_gene_short=[]
	true_y=pd.DataFrame()
	sm = ADASYN(random_state=42,n_neighbors=int(args.n_neighbors))

	X_sm, Y_sm = sm.fit_resample(X, Y)

	plt.figure(figsize=(8, 8))
	sns.countplot('Disease_Association', data=Y_sm)
	plt.title('Classes Distribution')
	plt.savefig('ClassesDistribution-aftersmote.pdf',bbox_inches='tight')
	plt.clf()
	for train_index, test_index in kf.split(X,Y):

		counter+=1
		X_train, X_validate = X.iloc[train_index], X.iloc[test_index]
		y_train, y_validate = Y.iloc[train_index], Y.iloc[test_index]
		print(y_train.Disease_Association.value_counts())
		X_train, y_train = sm.fit_resample(X_train, y_train)

		n_estimators = list(range(1,X_train.shape[1]))
		train_results=[]


		train_rf_predictions,train_rf_probs_0,train_rf_probs_1,rf_predictions,rf_probs_0,rf_probs_1 = cross_validate(X_train,X_validate,y_train,y_validate,counter,definitions)

		validate_gene_short.extend(list(X_validate.index.get_level_values(1)))
		validation.extend(rf_predictions)
		validation_prob_0.extend(rf_probs_0)
		validation_prob_1.extend(rf_probs_1)

		true_y=true_y.append(y_validate)


	test_rf_predictions,test_rf_probs_0,test_rf_probs_1 = predict(X_sm, Y_sm,X_test)

	df_test_pred = pd.DataFrame({'Gene': list(X_test.index.get_level_values(0)),'Gene_short_name' : list(X_test.index.get_level_values(1)),'predictions':test_rf_predictions, str(definitions[0]):test_rf_probs_0, str(definitions[1]):test_rf_probs_1})

	df_test_pred.to_csv('predictions.tsv',sep='\t')

	fpr, tpr, threshold =roc_curve(true_y["Disease_Association"].tolist(),validation_prob_1)
	i = np.arange(len(tpr))
	roc = pd.DataFrame({'fpr' : pd.Series(fpr, index=i),'tpr' : pd.Series(tpr, index = i), '1-fpr' : pd.Series(1-fpr, index = i), 'tf' : pd.Series(tpr - (1-fpr), index = i), 'thresholds' : pd.Series(threshold, index = i)})
	roc['gmean']=roc['tpr']*roc['1-fpr']
	ix = argmax(roc['gmean'])
	print('Threshold',roc.iloc[ix])
	roc.iloc[(roc.tf-0).abs().argsort()[:1]]
	roc.to_csv('roc.tsv', sep='\t')
	fig, ax = pl.subplots()
	pl.plot(roc['tpr'])
	pl.plot(roc['1-fpr'], color = 'red')
	pl.xlabel('1-False Positive Rate')
	pl.ylabel('True Positive Rate')
	pl.title('Receiver operating characteristic')
	ax.set_xticklabels([])
	pl.savefig('roc.png'+ str(counter)+'.pdf',bbox_inches='tight')

	df_valid = pd.DataFrame({'Gene_short_name' : validate_gene_short,'predictions': validation, str(definitions[0]):validation_prob_0, str(definitions[1]):validation_prob_1,'True':true_y["Disease_Association"].tolist()})


	df_valid.to_csv('predictions_validation.tsv', sep='\t')
	df_cm=pd.crosstab(df_valid['True'], df_valid['predictions'], rownames=['Actual Label'], colnames=['Predicted Label'])
	print(df_cm)
	df_valid_prob=df_valid

	candidategenes(df_test_pred,candidate_genes_df,candidate_genes)



def candidategenes(df_genes_pred,candidate_genes_df,candidate_genes):

	print(df_genes_pred['predictions'].value_counts())
	df_genes_pred.to_csv('All.tsv', sep='\t')
	df_final = candidate_genes_df.merge(df_genes_pred, on='Gene_short_name',how="left")
	print(df_final.loc[df_final['Gene_short_name'] == 'HMGB1','predictions'])
	print(df_final.loc[df_final['Gene_short_name'] == 'HOXD13','predictions'])
	print(df_final.loc[df_final['Gene_short_name'] == 'FBXW4','predictions'])
	print(df_final.loc[df_final['Gene_short_name'] == 'SHH','predictions'])

	df_final.to_csv(os.path.basename(candidate_genes), sep='\t')



def cross_validate(X_train,X_validate,y_train,y_validate,counter,definitions):

	model = RandomForestClassifier(n_estimators=int(args.n_estimators),  criterion='gini', max_depth=int(args.max_depth), min_samples_split=int(args.min_samples_split), min_samples_leaf=int(args.min_samples_leaf), min_weight_fraction_leaf=0.0, max_features=args.max_features, max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, bootstrap=bool(args.bootstrap), oob_score=True, n_jobs=None, random_state=None, verbose=0, warm_start=False, class_weight=None, ccp_alpha=0.0, max_samples=None)

	plt.clf()
	visualizer = ROCAUC(model)
	visualizer.fit(X_train, y_train)        # Fit the training data to the visualizer
	visualizer.score(X_validate, y_validate)        # Evaluate the model on the test data
	print('visualizer',visualizer.score(X_validate, y_validate))

	model.fit(X_train, y_train)
	n_nodes = []
	max_depths = []



	# Training predictions
	train_rf_predictions = model.predict(X_train)
	train_rf_probs_0 = model.predict_proba(X_train)[:, 0]
	train_rf_probs_1 = model.predict_proba(X_train)[:, 1]




	# Testing predictions
	rf_predictions = model.predict(X_validate)
	rf_probs_0 = model.predict_proba(X_validate)[:, 0]
	rf_probs_1 = model.predict_proba(X_validate)[:, 1]

	plt.clf()

	# Plot formatting
	reversefactor = dict(zip(range(3),definitions))
	y_validate_cm = np.vectorize(reversefactor.get)(y_validate).flatten()
	rf_predictions_cm = np.vectorize(reversefactor.get)(rf_predictions)

	df_cm=pd.crosstab(y_validate_cm, rf_predictions_cm, rownames=['Actual Label'], colnames=['Predicted Label'])
	df_cm.to_csv('confusion_matrix'+ str(counter)+'.tsv',sep='\t')

	plt.figure(figsize=(20,20))


	return train_rf_predictions,train_rf_probs_0,train_rf_probs_1,rf_predictions,rf_probs_0,rf_probs_1

def predict(X_train, y_train, X_test):
	model = RandomForestClassifier(n_estimators=int(args.n_estimators),  criterion='gini', max_depth=int(args.max_depth), min_samples_split=int(args.min_samples_split), min_samples_leaf=int(args.min_samples_leaf), min_weight_fraction_leaf=0.0, max_features=args.max_features, max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, bootstrap=bool(args.bootstrap), oob_score=True, n_jobs=None, random_state=None, verbose=0, warm_start=False, class_weight=None, ccp_alpha=0.0, max_samples=None)




	model.fit(X_train, y_train)
	importances = model.feature_importances_
	std = np.std([tree.feature_importances_ for tree in model.estimators_], axis=0)
	forest_importances = pd.Series(importances, index=X_train.columns)
	forest_importances.to_csv('forest_importances.tsv', sep='\t')
	fig, ax = plt.subplots()
	forest_importances.plot.bar(yerr=std, ax=ax)
	ax.set_title("Feature importances using MDI")
	ax.set_ylabel("Mean decrease in impurity")
	fig.tight_layout()
	plt.savefig('variableimp.pdf',bbox_inches='tight')
	test_rf_predictions = model.predict(X_test)
	test_rf_probs_0 = model.predict_proba(X_test)[:, 0]
	test_rf_probs_1 = model.predict_proba(X_test)[:, 1]

	return test_rf_predictions,test_rf_probs_0,test_rf_probs_1





if __name__ == "__main__":
	main()
