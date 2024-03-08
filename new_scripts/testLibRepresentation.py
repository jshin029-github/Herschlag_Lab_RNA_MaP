import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import seaborn as sns

bc = pd.read_csv("3chips-barcode-assignments.csv", sep='\t')
lib = pd.read_csv("rclib_unique.libchar", sep=",")

replib = lib.loc[ np.in1d(lib.loc[:, 'variant'].tolist(), bc.loc[:, 'variant'].tolist() ) ]

print "Found %d matched design sequences = %3.2f %%"%(replib.shape[0], 100*float(replib.shape[0])/lib.shape[0])

replib.to_csv("rclib_assigned.libchar", sep=",", header=True, index=False)

num_bcs =  bc.variant.value_counts()

num_bcs = pd.Series(num_bcs, name="No. of Assigned Barcodes")



