import sys
import math 
import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
from sklearn.linear_model import Lasso
from sklearn.linear_model import LinearRegression

from sklearn.model_selection import train_test_split


def main():
	feat=sys.argv[1]
	df = pd.read_csv(feat,sep='\t')
	#df = df.dropna()
	#df = df.loc[df['Disease-Association']]
	df = df.set_index(['hg_gene_stable_id','mm_gene_stable_id', 'Gene_symbol'])
	df = df[df['Disease-Association'].notna()]
	print(df)	#df=df.drop(['association_score.overall_y','association_score.datatypes.genetic_association','association_score.datatypes.somatic_mutation','association_score.datatypes.known_drug','association_score.datatypes.affected_pathway','association_score.datatypes.rna_expression','association_score.datatypes.literature','association_score.datatypes.animal_model','association_score.overall_x','HPOLabel','HPO_id','disease_ID'], axis=1)
	df = df.fillna(0)
	Y = df[['Disease-Association']]
	X = df.drop(['Disease-Association'],axis=1)
	X_train,X_test,y_train,y_test=train_test_split(X,Y, test_size=0.3, stratify=Y)
	lamdas = [0.000000001,0.00000001,0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1]
	lasso(X_train,X_test,y_train,y_test,lamdas)
	
def lasso(X_train,X_test,y_train,y_test,lamdas):
	l_num = len(lamdas)
	coefs=[]
	train_r_squared = np.zeros(l_num)
	test_r_squared = np.zeros(l_num)

	for ind, i in enumerate(lamdas):    
		reg = Lasso(alpha = i)
		reg.fit(X_train, y_train)
		train_r_squared[ind] = reg.score(X_train, y_train)
		test_r_squared[ind] = reg.score(X_test, y_test)
		coefs.append(reg.coef_)


	plt.figure(figsize=(18, 8))
	plt.plot(train_r_squared, 'bo-', label=r'$R^2$ Training set', color="darkblue", alpha=0.6, linewidth=3)
	plt.plot(test_r_squared, 'bo-', label=r'$R^2$ Validation set', color="darkred", alpha=0.6, linewidth=3)
	plt.xlabel('Lamda index'); plt.ylabel(r'$R^2$')
	plt.xlim(0, 8)
	plt.title(r'Evaluate lasso regression with lamdas: 0 = 0.000000001, 1 = 0.00000001, 2= 0.00000001, 3 = 0.000001, 4 = 0.00001, 5= 0.0001, 6= 0.001, 7 = 0.01, 8= 0.1, 9=1')
	plt.legend(loc='best')
	plt.grid()
	plt.savefig('lasso-R2-all3.pdf') 
	df_lam = pd.DataFrame(test_r_squared*100, columns=['R_squared'])
	df_lam['lambda'] = (lamdas)
	# returns the index of the row where column has maximum value.
	df_lam.loc[df_lam['R_squared'].idxmax()]
	test_r_squared = test_r_squared[::-1]
	print(test_r_squared,len(test_r_squared)-np.argmax(test_r_squared)-1)
	print('lamda ',lamdas[len(test_r_squared)-np.argmax(test_r_squared)-1])
	reg_best = Lasso(alpha = lamdas[np.argmax(test_r_squared)])
	reg_best.fit(X_train, y_train)
	reg_best.coef_
	best = pd.Series(reg_best.coef_, index=X_train.columns)
	best.to_csv('coefs.tsv', sep='\t')
	plt.figure(figsize=(15, 20))
	colormap = plt.cm.nipy_spectral
	colors = [colormap(i) for i in np.linspace(0, 1,len(X_test.columns))]
	ax = plt.gca()
	ax.set_prop_cycle('color', colors)
	ax.plot(lamdas, coefs)
	ax.set_xscale('log')
	ax.legend(X_test.columns, loc='lower right')
	plt.axis('tight')
	plt.xlabel('lambda')
	plt.ylabel('weights')
	plt.savefig('lasso-all3.pdf') 
	
if __name__ == "__main__":
    main()	
	
