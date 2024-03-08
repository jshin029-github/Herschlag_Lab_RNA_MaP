# Usage: python plotOverallFLuorQC.py [all.annot.green.CPseries] [all.annot.red.CPseries] [conditions.txt] [plotdir]
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import sys

if len(sys.argv) != 7:
    print "Usage python plotClusterHeatmap.py [input.CPseries] [input.Clus.gz] [Block X] [Block Y] [Filter Flag] [Plot Dir]"
    sys.exit(0)

# Read in inputs:
CPseries = pd.read_csv(gzip.open(sys.argv[1], 'rb'), sep='\t')
print CPseries.head(5)

Cluster = pd.read_csv(gzip.open(sys.argv[2], 'rb'), sep='\t')
print Cluster.head(5)

BlockX = int(sys.argv[3])
BlockY = int(sys.argv[4])
Filter = sys.argv[5]
PlotDir = sys.argv[6]

# Trim CPseries dataset
CPseries = CPseries.iloc[:, [0, -1]]
CPseries = CPseries.rename(cols={CPseries.columns[0]: 'clusterID', CPseries.columns[1]: 'Fluorescence'})
print CPseries.head(5)

# Filter Cluster dataset
Cluster = Cluster.loc[CLuster['clusterType'] == Filter]
X = Cluster['X']
Y = Cluster['Y']
F = Cluster['Fluorescence']

DF = pd.DataFrame()
DF['X'] = X
DF['Y'] = Y
DF['F'] = F
DF = DF.pivot(index='Y', columns="X", values="F")



# Join datasets:
df = Cluster.merge(CPseries, on = "clusterID", how="inner")
