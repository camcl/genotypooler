import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgba
import sys, os
import argparse
import numpy as np
import pandas as pd
from collections import Counter
from scipy.special import logit

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import dataframe as vcfdf

"""
Exploratory Data Analysis for the NIAB Diverse MAGIC wheat population

The data are:
- Genotypes for the 16 founders
- Genotypes for the inbred lines
- Genotypes are treated in GT format

Options for running the script:
- Paths to GT data (true or pooled, LD or LD + HD)
- Output directory for saving tables and figures
- Type of AF data to use (for binning): MAF or AAF

This scripts generates:
- SNPs counts per MAF or AAF bin as LaTeX tables (markers-bin-counts.tex)
- SNPs counts per MAF or AAF bin as jason files (markers-bin-counts.json)
- Heterozygote counts per MAF or AAF bin as LaTeX tables (het-bin-counts.tex)
- Heterozygote counts per MAF or AAF bin as jason files (het-bin-counts.json)
- counts of unassayed SNPs per MAF or AAF bin as LaTeX tables (missing-pooled-bin.tex)
  /!\ Only GT format should be used for computing missing rate!
- counts of unassayed SNPs per MAF or AAF bin as json files (missing-pooled-bin.json)
(- pooled MAF per true MAF or AAF bin as LaTeX tables (maf-pooled-bin.tex))
(- pooled MAF per true MAF or AAF bin as json files (maf-pooled-bin.json))

Command line usage (assuming the current directory is genotypooler/graphtools)
$ python3 -u eda_tables.py /home/camille/MagicWheat/data/20220611/1/Chr1.SNPs.pruned.nomiss.vcf.gz /home/camille/MagicWheat/src/genotypooler/examples/Chr1.SNPs.pruned.nomiss.pooled.imputed.vcf.gz /home/camille/MagicWheat/src/genotypooler/examples/results/beagle-no-map
"""

### COMMAND-LINE PARSING AND PARAMETERS
# TODO: GT or GL as cmd line param
parser = argparse.ArgumentParser(description='Exploratory Data Analysis with plots'
                                             'for the genetic structure in a population'
                                             'of simulated pooled genotypes')
parser.add_argument('truegenos', metavar='trueg', type=str, help='File with true genotypes GT', default=None)
parser.add_argument('othergenos', metavar='otherg', type=str, help='File with other genotypes GT or GL', default=None)
parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
parser.add_argument('mask', metavar='mask', type=str, help='File with filtered data (GP format)', default=None)
parser.add_argument('changed', metavar='chg', type=int, help='If True, mask unchanged priors', default=None)
argsin = parser.parse_args()

af_data = 'maf'

paths = {
    'true': argsin.truegenos,
    'other': argsin.othergenos
}

fmask = argsin.mask
changed = argsin.changed

# Create gbool mask for data

if fmask != '.':
    dfobj = vcfdf.PandasMixedVCF(fmask, format='GP', indextype='chrom:pos', mask=None)
    gdata = dfobj.genotypes()
    if changed:
        gbool = gdata.applymap(lambda x: True if x[0] is None else False).values  # ~gdata... && 'changed_priors' -> changed priors only
    else:
        gbool = ~gdata.applymap(lambda x: True if x[0] is None else False).values  # gdata... && 'changed_priors' -> unchanged priors only
    print('Number of visible values = ', gbool.sum())
else:
    gbool = None

outdir = argsin.outdir  # '/home/camille/MagicWheat/src/genotypooler/examples/results/beagle-no-map'
if not os.path.exists(outdir):
    os.makedirs(outdir)
print('\r\nFigures and tables will be saved in {}'.format(outdir).ljust(80, '.'))

# Plots' layout params

pts_sz = 10  # scatter dots size
pts_alpha = 1.0  # transparency
figsize = 8
labfont = 10

true_genos = [0.0, 1.0, 2.0]
true_labels = ['0/0', '0/1', '1/1']
pooled_genos = [0.0, 1.0, 2.0, -0.5, 0.5, -1.0]
pooled_labels = ['0/0', '0/1', '1/1', '0/.', './1', './.']
genocolors = ['#047495', '#00035b', '#748b97',  # full GT
              '#dbb40c', '#c65102', '#80013f'  # missing GT
              ]
genotricmap = ListedColormap([to_rgba(co) for co in genocolors[:3]], name='geno_tri_cmap')

x_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.94, 0.96, 0.98, 1.0]
lab_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85, 0.92, 0.95, 0.97, 0.99]
x2_bins = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]  # MAF
lab2_bins = [0.025, 0.075, 0.15, 0.25, 0.35, 0.45]  # MAF
if af_data == 'maf':
    x_bins = x2_bins
    lab_bins = lab2_bins
lab_fmt = ['{:.2f}-{:.2f}'.format(i, j) for i, j in zip(x_bins[:-1], x_bins[1:])]


# Read and process data

try:
    dftrue = vcfdf.PandasMixedVCF(paths['true'], format='GT')
except:
    dftrue = vcfdf.PandasMixedVCF(paths['true'], format='GL')
try:
    dfother = vcfdf.PandasMixedVCF(paths['other'], format='GT', mask=~gbool)
    print(dfother.genotypes().head())
except:
    dfother = vcfdf.PandasMixedVCF(paths['other'], format='GL', mask=~gbool)
    print(dfother.genotypes().head())

if af_data == 'aaf':
    df0 = dftrue.concatcols(
        [dftrue.af_info, dftrue.missing_rate, dftrue.aaf, dftrue.het_rate, dftrue.hom_alt_rate, dftrue.hom_ref_rate]
    )
    df1 = dfother.concatcols(
        [dftrue.af_info, dftrue.aaf.rename(columns={'aaf': 'true_af'}), dfother.missing_rate, dfother.aaf]
    )
    df2 = dfother.concatcols(
        [dfother.af_info, dftrue.aaf.rename(columns={'aaf': 'true_af'}), dfother.het_rate, dfother.hom_alt_rate,
         dfother.hom_ref_rate]
    )

if af_data == 'maf':
    df0 = dftrue.concatcols([dftrue.aaf, dftrue.missing_rate, dftrue.het_rate, dftrue.maf])  # dftrue.maf.rename(columns={'maf': 'true_af'}),
    df3 = dfother.concatcols([dfother.hom_alt_rate, dfother.hom_ref_rate])

print('\nTrue data:')
print(df0.head(10))
# Zoom-in on variants in a given MAF range
print(df0[0.02 < df0['aaf']][df0['aaf'] <= 0.05].iloc[:])
if df3 is not None:
    print('\nOther data:')
    print(df3.head(10))


# Basic statistics for true and masked data sets per bin

# Genotype counts (variants x samples) in true data

cnts = dict([(lab, 0) for lab in lab_fmt])  # lab_bins
mafbin = pd.cut(dftrue.maf.values.squeeze(), x_bins, labels=lab_fmt, include_lowest=True)  # lab_bins
cnts.update(Counter(mafbin.dropna()))
cnts = dict([(k, v * dftrue.genotypes().shape[1]) for k, v in cnts.items()])

# Table with true counts
dfbincnts = pd.DataFrame.from_records([cnts])
dfbincnts.index = ['true']
dfbincnts['Total'] = dftrue.genotypes().size  # .sum(axis=1)
dfbincnts.to_latex(buf=os.path.join(outdir, 'markers-bin-counts.tex'),
                   sparsify=True,
                   multirow=True,
                   caption='SNPs counts per {} bin'.format('MAF' if af_data == 'maf' else 'AAF'),
                   label='tab:markers-bin-counts')
dfbincnts.to_json(os.path.join(outdir, 'markers-bin-counts.json'),
                  orient='columns',
                  index=True
                  )
print('\n')
print(dfbincnts)

# Table with proportions (with respect to tot number and to the bin itself?)
dfbinprop = dfbincnts / dfbincnts['Total'][-1]
dfbinprop.to_latex(buf=os.path.join(outdir, 'markers-bin-prop.tex'),
                   sparsify=True,
                   multirow=True,
                   float_format="%.3f",
                   caption='SNPs proportions per {} bin with respect to the total number of SNPs on the genetic map'.format('MAF' if af_data == 'maf' else 'AAF'),
                   label='tab:markers-bin-prop')
dfbinprop.to_json(os.path.join(outdir, 'markers-bin-prop.json'),
                  orient='columns',
                  index=True
                  )

# Counts of visible genotypes (variants x samples) in masked data
truebin = pd.cut(df0.maf.values.squeeze(), x_bins, labels=lab_fmt, include_lowest=True)
df0['MAF-bin'] = truebin
df3['MAF-bin'] = truebin  # this assumes that the variants are the same in both data sets and sorted in the same order

macnts = dfother.genotypes().count(axis=1)
macnts.index.name = 'variants'
macnts.name = 'counts'
otherma = df3['MAF-bin'].to_frame()\
    .join(macnts, how='left')\
    .groupby(['MAF-bin'])['counts']\
    .sum()\
    .to_frame()\
    .transpose()
otherma.columns = otherma.columns.astype(str)
otherma['Total'] = otherma.loc['counts'].sum()

# Proportion of visible genotypes (variants x samples) in masked data w.r.t. the counts of genotypes per bin
maprops = otherma.divide(dfbincnts.values)\
    .round(decimals=3)\
    .rename({'counts': 'proportions'})

otherma = pd.concat([otherma, maprops], axis=0)
print('\n')
print(otherma)

otherma.to_latex(buf=os.path.join(outdir, 'unmasked-genotypes-bin.tex'),
                 sparsify=True,
                 multirow=True,
                 float_format="%.3f",
                 caption='Counts and proportions of unmasked genotypes per {} bin'.format('MAF' if af_data == 'maf' else 'AAF'),
                 label='tab:unmasked-genotypes-bin')
otherma.to_json(os.path.join(outdir, 'unmasked-genotypes-bin.json'),
                orient='columns',
                index=True
                )

