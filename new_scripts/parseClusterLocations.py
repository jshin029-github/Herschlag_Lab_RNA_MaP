import pandas as pd
import numpy as np
import gzip as gz
import os, sys

# Usage: python ... [input file] [output file]
#       input file: CPseq.gz file containing atleast cluster ID in first column
#       output file: output csv.gz file to write cluster ID and x , y columns (2d coordinates) and other info to

df = pd.read_csv(gz.open(sys.argv[1], 'rb'), header=None, sep="\t")
clusterIDs = df[df.columns[0]].tolist()
clusterTypes = df[df.columns[1]].tolist()

lane = []
tile = []
x = []
y = []

for clusterID in clusterIDs:
    words = clusterID.strip().split(':')
    L = int(words[3])
    T = int(words[4])
    X = int(words[5])
    Y = int(words[6])
    x.append(X)
    y.append(Y)
    lane.append(L)
    tile.append(T)

df = pd.DataFrame()
df['clusterID'] = clusterIDs
df['clusterType'] = clusterTypes
df['lane'] = lane
df['tile'] = tile
df['x'] = x
df['y'] = y

print df.head(10)
df.to_csv(sys.argv[2], sep="\t", header=True, index=False, compression="gzip")




