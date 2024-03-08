import pandas as pd
import numpy as np
import sys, gzip

# Usage: python filterBadClusters.py [Input.csv.gz] [BadClusters.csv] [Output.csv.gz]

inputfile = sys.argv[1]
filterfile = sys.argv[2]
outputfile = sys.argv[3]

df = pd.read_csv(gzip.open(inputfile, 'rb'), sep="\t")
print "Initially found: %d clusters"%(df.shape[0])
fl = pd.read_csv(filterfile, sep=",")

clusID = df.columns[0]

badclusters = fl['clusterID']
dfl = df[~df[clusID].isin(badclusters)]
print "Finally have : %d clusters"% (dfl.shape[0])
dfl.to_csv(outputfile, sep="\t", compression="gzip", index=False)
