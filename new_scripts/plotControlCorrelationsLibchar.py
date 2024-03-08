# Python script to compare data and plot with values in ref_data via 
# Plots generated in ./plots and named: plots/ControlCorr_*.png
# Usage: python plotControlCorrelations.py Exp.CPvariant.gz Ref.CPvariant VarianMap.csv Suffix

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import gzip
import sys
import numpy as np

RT = 0.582
conc = 1E-9

if len(sys.argv) != 5:
    print 'Usage: python plotControlCorrelations.py Exp.CPvariant.gz Ref.CPvariant VarianMap.csv Suffix'
    sys.exit(0)

# Read Exp CPvariant file - file of bootstrapped fits, and format it
print 'Reading in Exp data'
df = pd.read_csv(gzip.open(sys.argv[1], "rb"), sep="\t")
df = df.rename(columns = {'Unnamed: 0': 'variant'})

# Read Ref CPvariant file - file of bootstrapped fits, and format it
print 'Reading in Ref data'
refdf = pd.read_csv(gzip.open(sys.argv[2], "rb"), sep="\t")
refdf = refdf.rename(columns = {'Unnamed: 0': 'ref_variant', 'dG':'dG_ref'})

# Read variant map csv file
print "Reading in Variant Map"
varmap = pd.read_csv(sys.argv[3], sep="\t")
varmap = varmap.loc[:, ['ref_variant', 'exp_variant']]

# Merge everything together to plot in df
print "Merging data together to plot"
df = df.merge(varmap, left_on = "variant", right_on="exp_variant", how="inner")
df = df.merge(refdf, left_on = "ref_variant", right_on="ref_variant", how="inner")

# Plot dG Exp vs dG Reference:
print "Plotting dG correlation"
sns.lmplot(data=df, x='dG_ref', y='dG', markers='o',scatter_kws={'s': 20})
R2 = np.corrcoef(df['dG_ref'], df['dG'])[0,1]
print "     Plotting dG Exp vs dG Ref"
print "     Found %d points"%(df.shape[0])
print "     Obtained R^2 = ", R2
print "     Obtained RMSE = ", np.sqrt(np.mean((df['dG_ref']-df['dG'])**2))
plt.xlabel('Reference dG')
plt.ylabel('Measured dG')
plt.savefig('plots/ControlCorr_dG_correlation_%s.png'%(sys.argv[-1]))
plt.clf()
    

# Plot logKd Exp vs logKd Reference (Kd < 10^8)
print "Plotting Kd correlation"
df['logKd'] = np.log10(np.exp(df['dG']/RT)/conc)
df['logKd_ref'] = np.log10(np.exp(df['dG_ref']/RT)/conc)
sns.lmplot(data=df, y='logKd', x='logKd_ref', fit_reg=False, markers="o", scatter_kws={'s':20})
plt.xlim([0,6])
plt.ylim([0,6])
plt.xlabel('Reference log10(Kd)')
plt.ylabel('Measured log10(Kd)')
plt.savefig('plots/ControlCorr_Kd_correlation_%s.png'%(sys.argv[-1]))
plt.clf()

