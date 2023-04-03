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

# TODO: Separate script for plotting aaf-diff and script for LaTeX tables
# TODO: Edit descriptions of scripts
# TODO: bootrapped T-test for significance? (maybe overkilled?)
# TODO: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html, https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html


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

Command line usage (assuming the current directory is genotypooler/sandbox)
$ python3 -u eda_tables.py /home/camille/MagicWheat/data/20220611/1/Chr1.SNPs.pruned.nomiss.vcf.gz /home/camille/MagicWheat/src/genotypooler/examples/Chr1.SNPs.pruned.nomiss.pooled.imputed.vcf.gz /home/camille/MagicWheat/src/genotypooler/examples/Chr1.SNPs.pruned.nomiss.pooled.imputed.vcf.gz /home/camille/MagicWheat/src/genotypooler/examples/results/beagle-no-map
"""

### COMMAND-LINE PARSING AND PARAMETERS
# TODO: GT or GL as cmd line param
parser = argparse.ArgumentParser(description='Exploratory Data Analysis with plots'
                                             'for the genetic structure in a population'
                                             'of simulated pooled genotypes')
parser.add_argument('truegenos', metavar='trueg', type=str, help='File with true genotypes GT', default=None)
parser.add_argument('othergenos', metavar='otherg', type=str, help='File with other genotypes GT or GL', default=None)
parser.add_argument('outdir', metavar='outdir', type=str, help='Directory to save the plots', default=None)
argsin = parser.parse_args()

af_data = 'maf'

paths = {
    'true': argsin.truegenos,
    'other': argsin.othergenos
}
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
    dfother = vcfdf.PandasMixedVCF(paths['other'], format='GT')
    print(dfother.genotypes().head())
except:
    dfother = vcfdf.PandasMixedVCF(paths['other'], format='GL')
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
    df0 = dftrue.concatcols([dftrue.aaf, dftrue.maf.rename(columns={'maf': 'true_af'}), dftrue.missing_rate, dftrue.het_rate, dftrue.maf])
    df3 = dfother.concatcols([dfother.aaf, dfother.maf, dfother.missing_rate, dfother.het_rate, dfother.hom_alt_rate, dfother.hom_ref_rate])

print('\nTrue data:')
print(df0.head(10))
if df3 is not None:
    print('\nOther data:')
    print(df3.head(10))

print('\r\nOverall missing rates:')
tottruemiss = dftrue.missing_rate.mean() * 100
totothermiss = dfother.missing_rate.mean() * 100
print('True: ', tottruemiss)
print('Other: ', totothermiss)

# Basic statistics for LD and HD data sets per bin

# Markers counts

cnts = dict([(lab, 0) for lab in lab_fmt])  # lab_bins
hdbin = pd.cut(dftrue.maf.values.squeeze(), x_bins, labels=lab_fmt, include_lowest=True)  # lab_bins
cnts.update(Counter(hdbin.dropna()))
hdcnts = cnts

# Table with counts
dfbincnts = pd.DataFrame.from_records([hdcnts])
dfbincnts.index = ['true']
dfbincnts['Total'] = dfbincnts.sum(axis=1)
dfbincnts.to_latex(buf=os.path.join(outdir, 'markers-bin-counts.tex'),
                   sparsify=True,
                   multirow=True,
                   caption='SNPs counts per {} bin'.format('MAF' if af_data == 'maf' else 'AAF'),
                   label='tab:markers-bin-counts')
dfbincnts.to_json(os.path.join(outdir, 'markers-bin-counts.json'),
                  orient='columns',
                  index=True
                  )

# Table with proportions
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

# Assign variants to bins

truebin = pd.cut(df0.true_af.values.squeeze(), x_bins, labels=lab_fmt, include_lowest=True)
df0['AF-bin'] = truebin
df3['AF-bin'] = truebin  # this assumes that the variants are the same in both data sets and sorted in the same order

# Missing counts

truemiss = df0.groupby(['AF-bin'])['missing_rate'].mean().to_frame() * 100  # display as percentage
othermiss = df3.groupby(['AF-bin'])['missing_rate'].mean().to_frame() * 100

binmiss = truemiss.join(othermiss, how='left', lsuffix='_true', rsuffix='_other').transpose()
binmiss.columns = binmiss.columns.to_list()
binmiss['Total'] = np.concatenate([tottruemiss.values, totothermiss.values])
binmiss.index = ['true', 'other']
print('\r\nMissing rate per bin:')
print(binmiss)

binmiss.to_latex(buf=os.path.join(outdir, 'missing-other-bin.tex'),
                 sparsify=True,
                 multirow=True,
                 float_format="%.3f",
                 caption='Proportion of missing genotypes per {} bin'.format('MAF' if af_data == 'maf' else 'AAF'),
                 label='tab:missing-other-bin')
binmiss.to_json(os.path.join(outdir, 'missing-other-bin.json'),
                orient='columns',
                index=True
                )


# Heterozygote counts

truehet = df0.groupby(['AF-bin'])['het_rate'].mean().to_frame() * 100
otherhet = df3.groupby(['AF-bin'])['het_rate'].mean().to_frame() * 100

binhet = truemiss.join(otherhet, how='left', lsuffix='_true', rsuffix='_other').transpose()
binhet.columns = binhet.columns.to_list()
binhet['Total'] = np.concatenate([truehet.sum().values, otherhet.sum().values])
binhet.index = ['true', 'other']
print('\r\nHeterozygosity rate per bin:')
print(binhet)

binhet.to_latex(buf=os.path.join(outdir, 'het-other-bin.tex'),
                sparsify=True,
                multirow=True,
                float_format="%.3f",
                caption='Proportion of heterozygous genotypes per {} bin'.format('MAF' if af_data == 'maf' else 'AAF'),
                label='tab:het-other-bin')
binhet.to_json(os.path.join(outdir, 'het-other-bin.json'),
               orient='columns',
               index=True
               )
