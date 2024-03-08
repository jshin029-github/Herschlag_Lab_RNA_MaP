# Python script too plot and save info and stats on the different 
# Libraries present in Kds fitted
# Plots generated in ./plots and named: plots/LibStats_*.png
# Usage: python plotLibraryStats.py Bootstrap.CPvariant.gz Unique.libChar outputdir

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import sys
import numpy as np
import os

RT = 0.582
conc = 1E-9

if len(sys.argv) != 4:
    print "Usage: python plotLibraryStats.py Bootstrap.CPvariant.gz Unique.libChar outputdir"
    sys.exit(0)

# Read CPvariant file - file of bootstrapped fits, and format it
print 'Reading in CPvariant data'
df = pd.read_csv(gzip.open(sys.argv[1], "rb"), sep="\t")
df = df.rename(columns = {'Unnamed: 0': 'variant'})

# Read LibChar file and format it
print 'Reading in LibChar data'
lc = pd.read_csv(sys.argv[2])
#### NOTE: Change this to group by different columns in your libchar!
project_name = 'body_plan'
lc_colnames = lc.columns
final_colnames = []
for name in lc_colnames:
    if ('Unnamed' not in name) and (name != 'rcseq'):
        final_colnames.append(name)
lc = lc.loc[:, final_colnames]
lc = lc.rename(columns = {'variant':'libvariant'})

# Combine CPvariant with LibChar
print 'Combining LibChar with CPvariant file'
df = df.merge(lc, left_on='variant', right_on='libvariant', how="left")
del df['libvariant']

# Split up into sublibraries
print 'Finding projects present in CPvariant data'
splitdf = {}
sublibs = df[project_name].unique()
print 'Found %d unique projects / sublibaries'%(sublibs.shape[0])

print 'Generating sublibrary split results and stats'
count = 0
fout = open(os.path.join(sys.argv[3], 'LibraryStats.txt'), "w")

fout.write("%30s\t%6s\t%6s\t%4s\t%6s\t%4s\t%4s\n"%("LibraryName", "NVar", "TotVar", "%Lib", "NC>5", "MedNC", "MeddG"))

for sublib in sublibs:
    count += 1
    print "     Generating sublibrary %d / %d : %s "%(count, sublibs.shape[0], sublib)
    splitdf = df.loc[df[project_name] == sublib]
    nvars = splitdf.shape[0]
    print "         found %d variants"%(nvars)
    totalvars = lc.loc[lc[project_name] == sublib].shape[0]
    ncluster5 = np.sum(splitdf['num_tests'] > 5)
    nclustermedian = splitdf['num_tests'].median()
    quantiles = splitdf['dG'].quantile(q=[0.5])
    percent = float(nvars)*100.0/float(totalvars)
    fout.write("%30s\t%6d\t%6d\t%4.1f\t%6d\t%4.1f\t%4.2f\n"%(sublib, nvars, totalvars, percent, ncluster5, nclustermedian, quantiles[0.5]))
    outfilename = os.path.join(sys.argv[3], "LibSplit_%d.CPvariant"%(count))
    splitdf.to_csv(outfilename, sep="\t", index=False, header=True)












