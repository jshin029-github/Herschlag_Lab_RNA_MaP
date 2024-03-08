import pandas as pd, gzip
import sys

inputfile = sys.argv[1]
colname = sys.argv[2]
outputfile = sys.argv[3]

df = pd.read_csv(gzip.open(inputfile, 'rb'), sep='\t')
print df.shape
dff= pd.DataFrame()
dff['clusterID'] = df['clusterID']

for i in range(8):
    dff[str(i)] = df[colname]

print dff.head(5)
print dff.shape
dff.to_csv(outputfile, sep='\t', index=False, compression="gzip")
