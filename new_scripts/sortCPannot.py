import pandas as pd

cpannot = pd.read_csv("3chips_anyRNA.CPannot", sep='\t')
cpannot_sorted = cpannot.sort_values(by = ['variant_number', 'clusterID'])
cpannot_sorted.to_csv("3chips_anyRNA_sorted.CPannot", sep='\t', header=True, index=False)

