# Usage: python annotateClusters.py [CPseq.gz file] [barcode_assignment.csv] [output CPannot] [barcode col]
import numpy as np, pandas as pd, sys, gzip as gz

if len(sys.argv) != 5:
    print "python annotateClusters.py [CPseq.gz file] [barcode_assignment.csv] [output CPannot] [barcode col]:"
    sys.exit(0)

if sys.argv[1][-2:] == 'gz':
    cpseq = pd.read_csv(gz.open(sys.argv[1], 'rb'), sep="\t")
elif sys.argv[1][-5:] == 'CPseq':
    cpseq = pd.read_csv(sys.argv[1], sep="\t")
else:
    sys.exit('please give a .CPseq or .CPseq.gz file')
outfile = sys.argv[-2]
barcodecol = int(sys.argv[-1])-1

cpseqreduced = cpseq.iloc[:, [0, barcodecol]]


cpseqreduced.columns = ['clusterID', 'barcode']
cpseqreduced = cpseqreduced.sort_values("barcode")
bc = pd.read_csv(sys.argv[2], sep="\t")
bcreduced = bc.iloc[:, [0,1]]

print "Read in CPseq file: "
print "Total number of clusters sequenced : %d" %(cpseqreduced.shape[0])

print "Debug: Reduced CPseq DF"
print cpseqreduced.head()

print "Read in Barcode Assignment file: "
print "Total number of barcode assignments : %d" %(bcreduced.shape[0])

print "Debug: Reduced BCseq DF"
print bcreduced.head()

print "Assigning variants to clusters :"
df = cpseqreduced.merge(bcreduced, on="barcode", how="inner") 
df = df.sort_values(by="variant")
print "No. of clusters annotated: %d "%(df.shape[0])

print "Debug: Assigned cluster DF"
del df['barcode']
print df.head()
df = df.rename(columns={'variant':'variant_number'})
df.to_csv(outfile, sep="\t", compression="gzip", index=False, header=True)
