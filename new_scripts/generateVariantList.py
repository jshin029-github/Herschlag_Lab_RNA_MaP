# Usage: python generateVariantList.py [CPvariant] [outdir]
# Generates a list of Variant IDs (split into 10 files)

import numpy as np, pandas as pd, os, sys
import gzip

if len(sys.argv) != 3:
    print "Usage: python generateVariantList.py [Cpavariant] [outdir]"
    sys.exit(0)

cpvar = sys.argv[1]
outdir = sys.argv[2]

df = pd.read_csv(gzip.open(cpvar, 'rb'), sep='\t')
df = df.rename(columns={'Unnamed: 0': 'variant'})
variantIDs = np.array(df['variant'], dtype=int)

print "Found %d variants in CPvariant file"%(variantIDs.shape[0])
print "Splitting into variant ID files:"
array_list = np.split(variantIDs, 10)
counter = 1

for array in array_list:
    filename = "%s/variant_list_%d.txt"%(outdir, counter)
    print "\tWriting variant IDs: %d ... %d to %s"%(array[0], array[-1], filename)
    np.savetxt(filename, array, newline="\t", fmt="%6d")
    counter += 1
