import pandas as pd
df = pd.read_csv('Projects/Prioritization/monarch/gene_phenotype.all.tsv', sep='\t') #downloaded from monarch

df1 = df[df.subject_taxon_label.isin(['Homo sapiens','Mus musculus'])]

#df1.to_csv('gene_phenotype.human_mouse.tsv', sep='\t')


df2= pd.read_csv('Projects/Prioritization/monarch/Cardiac_MPhenotype_MGenotype.tsv', sep='\t') # downloaded from mousemine https://www.mousemine.org/mousemine/begin.do


cardiac_phenotype = df2.Ontology_Name.to_list()

df3=df1[df1.object_label.isin(cardiac_phenotype)]

df3 = df3.groupby(['subject_label','subject_taxon_label']).agg(lambda x:x.tolist())
df3.to_csv('Projects/Prioritization/monarch/monarch_cardiac_genes.tsv', sep='\t')

df4= pd.read_csv('Projects/Prioritization/monarch/Limb_MPhenotype_MGenotype.tsv', sep='\t') # downloaded from mousemine https://www.mousemine.org/mousemine/begin.do


limb_phenotype = df4.Ontology_Name.to_list()

df5=df1[df1.object_label.isin(limb_phenotype)]

df5 = df5.groupby(['subject_label','subject_taxon_label']).agg(lambda x:x.tolist())
df5.to_csv('Projects/Prioritization/monarch/monarch_cardiac_genes.tsv', sep='\t')
