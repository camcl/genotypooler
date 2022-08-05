"""
How does pooling decode the data? How "bad" does pooling make the data? What kind of missing data given the allele frequency?

Usage example:
$ python3 -u genotypes_decoding.py /home/camille/PoolImpHuman/data/20200812/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.pooled.snps.gt.vcf.gz /home/camille/PoolImpHuman/results/20200812
$ python3 -u genotypes_decoding.py /home/camille/PoolImpHuman/data/20200812/IMP.chr20.snps.gt.vcf.gz /home/camille/PoolImpHuman/data/20200812/IMP.chr20.missingHD.fullLD.snps.gt.vcf.gz /home/camille/PoolImpHuman/results/20200812
"""

import numpy as np
import pandas as pd
from sklearn.metrics import multilabel_confusion_matrix, confusion_matrix
from sklearn.preprocessing import *
from scipy.stats import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgba_array, to_rgba
import seaborn as sns
import os, sys
import argparse

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import dataframe as vcfdf


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Barplots for genotypes states in pooled data')
parser.add_argument('truefile', metavar='truef', type=str, help='File with true genotypes GT', default=None)
parser.add_argument('pooledfile', metavar='pooledf', type=str, help='File with pooled genotypes decoded into GT', default=None)
parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
argsin = parser.parse_args()

truef = argsin.truefile
# argsin.truefile
# '/home/camille/1000Genomes/src/genotypooler/examples/IMP.chr20.snps.gt.vcf.gz'
pooledf = argsin.pooledfile
# argsin.pooledfile
# '/home/camille/1000Genomes/src/genotypooler/examples/IMP.chr20.pooled.snps.gt.vcf.gz'
outdir = argsin.outdir
if not os.path.exists(outdir):
    os.mkdir(outdir)
print('\r\nFigures will be saved in {}'.format(outdir).ljust(80, '.'))


# Plotting features and parameters
true_genos = [0.0, 1.0, 2.0]
true_labels = ['0/0', '0/1', '1/1']
pooled_genos = [0.0, 1.0, 2.0, -0.5, 0.5, -1.0]
pooled_labels = ['0/0', '0/1', '1/1', '0/.', './1', './.']
mMlabels = ['M/M', 'M/m', 'm/m', 'M/.', './m', './.']  # minor-major allele

# dashes_styles = [(0, ()), (0, ()), (0, ()),  # full GT
#                  (0, (5, 5)), (0, (5, 5)), (0, (1, 1))  # missing GT
#                  ]
dashes_styles = ['-', '-', '-',  # full GT
                 '--', '--', '-.'  # missing GT
                 ]
barcolors = ['#047495', '#00035b', '#748b97',  # full GT
             '#dbb40c', '#c65102', '#80013f'  # missing GT
             ]
barcmap = ListedColormap([to_rgba(co) for co in barcolors])

xtype = 'maf'
#bin_step = 0.04
# x_bins = np.arange(0.0, 1.0 + bin_step, bin_step).round(decimals=2)
x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
x2_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.5]  # MAF
# lab_bins = np.arange(bin_step/2, 1.0, bin_step).round(decimals=2)
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]
lab_text = ['-'.join([str(x[0]), str(x[1])]) for x in list(zip(x_bins[:-1], x_bins[1:]))]
lab2_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.45]  # MAF
lab2_text = ['-'.join([str(x[0]), str(x[1])]) for x in list(zip(x2_bins[:-1], x2_bins[1:]))]

figsize=4
plt.rcParams["figure.figsize"] = [figsize*3, figsize + 1] if xtype == 'aaf' else [figsize*3, figsize * 2]

titlesz = 24
axlabsz= 20
axticksz = 16
legsz = 20
sns.set(font_scale=1.75)  # multiplication factor!

print('\r\nCounting genotypes'.ljust(80, '.'))


# Read and process data
print('\r\nReading data from {} and {}'.format(truef, pooledf).ljust(80, '.'))
dftrue = vcfdf.PandasMixedVCF(truef, format='GT')

try:
    dfpooled = vcfdf.PandasMixedVCF(pooledf, format='GL')
    pooled = dfpooled.trinary_encoding()
except KeyError:
    dfpooled = vcfdf.PandasMixedVCF(pooledf, format='GT')
    pooled = dfpooled.trinary_encoding()

n_markers, n_samples = dftrue.genotypes().shape


if xtype == 'aaf':
    af_bins = pd.cut(dftrue.aaf.values.squeeze(), bins=x_bins, labels=lab_bins, include_lowest=True)
    binned_af = pd.Series(af_bins, index=dftrue.variants, name='binned_aaf')

    pooled = pooled.join(binned_af)
    pooled['dataset'] = ['pooled'] * pooled.shape[0]
    pooled = pooled.reset_index().set_index(['variants', 'binned_aaf', 'dataset'])


if xtype == 'maf':
    maf_bins = pd.cut(dftrue.maf.values.squeeze(), bins=x2_bins, labels=lab2_bins, include_lowest=True)
    binned_maf = pd.Series(maf_bins, index=dftrue.variants, name='binned_maf')

    pooled = pooled.join(binned_maf)
    pooled['dataset'] = ['pooled'] * pooled.shape[0]
    pooled = pooled.reset_index().set_index(['variants', 'binned_maf', 'dataset'])

# Initialize counts for each AF bin and each genotype
print('\r\nCounting genotypes'.ljust(80, '.'))

if xtype == 'aaf':
    dfcounts = pd.DataFrame(data=None, index=lab_bins, columns=pooled_genos)

if xtype == 'maf':
    dfcounts = pd.DataFrame(data=None, index=lab2_bins, columns=pooled_genos)

for i in dfcounts.index:
    for j in dfcounts.columns:
        dfbins = pooled.loc[pooled.index.get_level_values('binned_{}'.format(xtype)) == i]
        dfbins.reset_index(inplace=True, drop=True)
        counts_geno = dfbins.where(dfbins == j, axis=0).count()
        dfcounts.loc[i, j] = counts_geno.sum()

if xtype == 'aaf':  # change x-ticks labels
    dfcounts.index = lab_text
    dfcounts.columns = pooled_labels

if xtype == 'maf':  # change x-ticks labels
    dfcounts.index = lab2_text
    dfcounts.columns = mMlabels


# scale result per bin
binscales = dfcounts.sum(axis=1)
dfcounts_scaled = pd.DataFrame(data=dfcounts.div(binscales, axis=0),  # minmax_scale(dfcounts.values, axis=0),
                               index=dfcounts.index,
                               columns=dfcounts.columns)
print(dfcounts_scaled)

dfcounts_sized = dfcounts / (n_samples * n_markers)


# Plot processed data
print('\r\nPlotting results'.ljust(80, '.'))
ax = dfcounts_scaled.plot(kind='bar', stacked=True, rot=45,
                          color=barcolors, style=dashes_styles)  # cmap = sns.set_palette('GnBu_d')
ax.set_xlabel('True {} allele frequency'.format('alternate' if xtype == 'aaf' else 'minor'), fontsize=axlabsz)
ax.set_ylabel('Proportions of genotypes scaled per {}-bin'.format(xtype.upper()), fontsize=axlabsz)
plt.title('Genotypes proportions in the study population', fontsize=titlesz)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'genotypes_hexa_scaled_proportions.pdf'))
plt.show()

ax_scaled = dfcounts_sized.plot(kind='bar', stacked=True, rot=45,
                                color=barcolors, style=dashes_styles)
ax.set_xlabel('True {} allele frequency'.format('alternate' if xtype == 'aaf' else 'minor'), fontsize=axlabsz)
ax_scaled.set_ylabel('Proportion of genotypes', fontsize=axlabsz)
plt.title('Genotypes proportions in the study population (total number of genotypes = {})'.format(n_samples * n_markers),
          fontsize=titlesz)
plt.tight_layout()
plt.savefig(os.path.join(outdir, 'genotypes_hexa_proportions.pdf'))
plt.show()
