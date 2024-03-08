# Combines a reference data set with your libChar file to give  a new formatted 
# dataset called CPref
# Requires: input dataset have a column called 'sequence' for comparison and it
#	    retains all fields with 'dG' in it from your reference dataset
# Usage: python processReferenceLibrary.py my.libChar ref.csv out.csv
# NOTE: Look throught comments for changes to make for ref lichar

import pandas as pd
import numpy as np
import sys

if len(sys.argv) != 4:
    print "Usage: python processReferenceLibrary.py my.libChar ref.csv out.csv"
    sys.exit(0)

lc = pd.read_csv(sys.argv[1])
df = pd.read_csv(sys.argv[2])
outfile = sys.argv[3]

# Filter out cols in libChar file and drop duplicates
lc = lc.loc[:, ['variant', 'chip_sequence']]
lc = lc.drop_duplicates('chip_sequence')
lc = lc.rename(columns = {'variant': 'exp_variant'})

# Filter out cols in ref dataframe and drop duplicates
dfColsToRetain = ['sequence']

# Uncomment for use with refLibchar
# dfColsToRetain = ['sequence', 'variant']

for colname in df.columns:
    if 'dG' in colname:
        dfColsToRetain.append(colname)
df = df.loc[:, dfColsToRetain]
df = df.drop_duplicates('sequence')

# Uncomment if using reference libchar
# df = df.rename(columns = {'variant': 'ref_variant'}) 

# Mask out common sequences and merge and write output
cmask = np.in1d(lc['chip_sequence'], df['sequence'])
lc_common = lc.loc[cmask]
df = df.merge(lc_common, how="inner", left_on="sequence", right_on="chip_sequence")
print "Found %d common sequences"%(df.shape[0])
del df['chip_sequence']
df.to_csv(outfile, sep='\t', index=False, header=True)



