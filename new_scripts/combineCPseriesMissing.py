## Usage: python combineCPseries.py [N tiles] [N concs] [all CPseries Dirs] [imagename template] [AllOutFilename] [anyRNAOutFilename] [CPannot.gz]
import sys, os
import numpy as np
import pandas as pd
import gzip

if len(sys.argv) < 7:
    print "Usage: python combineCPseries.py [N tiles] [N concs] [all CPseries Dirs] [imagename template] [AllOutFilename] [anyRNAOutFilename] [CPannot.gz]"
    sys.exit(0)

Ntiles = int(sys.argv[1])
imagetemp = sys.argv[-4]

TileFilenames = []
for i in range(1, Ntiles+1):
    #TileFilenames.append("fid_anyRNA_JL4CY_tile%03d_Bottom_filtered.CPseries"%(i))
    TileFilenames.append(imagetemp%(i))
print "Tiled Filenames:"
for n in TileFilenames:
    print "\t%s"%(n)

Nconcs = int(sys.argv[2])

CPseriesDirs = []
for j in range(Nconcs):
    CPseriesDirs.append(sys.argv[3+j])
print "CPseries Dir Filenames:"
for n in CPseriesDirs:
    print "\t%s"%(n)

print "\n\n"

outDir = sys.argv[-3]
anyRNAoutDir = sys.argv[-2]
cpannot = pd.read_csv(gzip.open(sys.argv[-1], 'rb'), sep="\t")
print cpannot.head(10)

# Concat all data frames vertically first:
DFseries = []
cpcount = 0

print " ***** Merging Vertically ***** "

for CPseriesDir in CPseriesDirs:
    print "Reading tiles from : %s" % (CPseriesDir)
    cpcount += 1
    # For each conc. point => read all available tile CPseries files and stack vertically
    firstPass = True
    DF  = None
    for TileFilename in TileFilenames:
	# If file doens't exist - continue on to next file
	    if not os.path.exists(os.path.join(CPseriesDir, TileFilename)):
	        print  "	File not found : %s " % (TileFilename)
   	        print "	Skipping processing steps"

	    else:
	        # Otherwise, files exists => Read and format it;
	        print "		Reading in file: %s "%(TileFilename)
	        d = pd.read_csv(os.path.join(CPseriesDir, TileFilename), sep='\t', header=None)
	        #d = d.loc[d[1] == 'anyRNA']
	        d = d.iloc[:, [0,1, -1]]
	        d.columns = ['clusterID','clusterFlag' , cpcount-1]
	        # File exists and is the first file in line
	        if  firstPass: ## Fist tile => read into DF
	            firstPass = False
	            #d = d.iloc[:, [0,1, -1]]
	            #d.columns = ['clusterID','clusterFlag' , cpcount-1]
	            DF = pd.DataFrame(d)
	        else:
	            #d = d.iloc[:, [0, -1]]
	            #d.columns = ['clusterID' , cpcount-1]
	            DF = pd.concat([DF, d])

    try:
        #DF.set_index('clusterID')
        print "		Final nrows = %d\n\n"%(DF.shape[0])
        print DF.head(20)
        print " -------------------------------------------- "
        DFseries.append(DF)
    except AttributeError:
        print "%s has no tiles" % CPseriesDir

print "No of dfs: ", len(DFseries)

# Merge vertically
print " ***** Merging Horizontally ***** "
DFCombined = pd.DataFrame(DFseries[0])

for i in range(1, len(DFseries)):
    print "Merging conc point %d"%(i)
    DFseries[i] = DFseries[i].iloc[:, [0, -1]]
    DFCombined = pd.merge(DFCombined, DFseries[i], how='outer', on='clusterID')

# Generated Full CPseries
print "Final Merged CPseries file: "
print DFCombined.head(20)

print "\n\n"

print "CPannot Data: "
print cpannot.head(20)

print "\n\n"

print "Combining CPannot data into CPseries files (assigning variant IDs):"
DFCombined = pd.merge(DFCombined, cpannot, how='left', on="clusterID")

print "\n\n"

print "Final combined and annotated CPseries file:"
print DFCombined.head(20)

print "\n\nGenerating all cluster output file %s"%(outDir)
DFCombined.to_csv(outDir, sep='\t', header=True, index=False, compression="gzip")

print "\n\nGenerating anyRNA cluster output file %s"%(anyRNAoutDir)
DFCombined = DFCombined.loc[DFCombined['clusterFlag'] == 'anyRNA']
del DFCombined['clusterFlag']
del DFCombined['variant_number']
DFCombined.to_csv(anyRNAoutDir, sep='\t', header=True, index=False, compression="gzip")
