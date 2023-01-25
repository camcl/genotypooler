import os, sys
import argparse
import numpy as np
import pandas as pd
import timeit

from typing import *

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs.metrics import quality as qual
from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.persotools.files import *

'''
This script addresses the following review for the manuscript:
Present in the same table: 
1. How many marker genotypes are know before imputation in both scenarios? (fully determined after pooling, or genotyped on LD map for the classical case),
2. How many markers are correctly imputed (exact matches), both with Phaser and with Beagle?

Arguments are parsed from a file where they must be written in this order:
pathout <>
truegt <>
truegl <>
pooledgl <>
pooledgt <<leave_blank_line_if_none>
imp1 <path/to/file/imputed/with/Beagle>
imp2 <path/to/file/imputed/with/Phaser>
date <>
bins <>
compute <>

Beware that strings for paths should be written just as text (without quotes!) in the argsfile!

Command line usage (assuming the current directory is genotypooler/graphtools)
$ python3 -u latex_tables_decoded_imputed.py @tablesargsfile.txt
'''

### COMMAND-LINE PARSING AND PARAMETERS (arguments are parsed from a file
parser = argparse.ArgumentParser(description='Plots the imputation performance'
                                             '(concordance and cross-entropy)',
                                 fromfile_prefix_chars='@')
parser.add_argument('pathout', metavar='out', type=str, help='Results directory', default=None)
parser.add_argument('truegt', metavar='gt', type=str, help='File with true data (GT format)', default=None)
parser.add_argument('truegl', metavar='gl', type=str, help='File with true data (GL format)', default=None)
parser.add_argument('pooledgl', metavar='poogl', type=str, help='File with pooled data (GL format)', default=None)
parser.add_argument('imp1', metavar='imp1', type=str, help='File with imputed data 1 (GT and/or GP formats)', default=None)
parser.add_argument('imp2', metavar='imp2', type=str, help='File with imputed data 2 (GT and/or GP formats)', default=None)
parser.add_argument('date', metavar='date', type=str, help='Date of the experiment (YYYMMDD)', default=None)
parser.add_argument('bins', metavar='bin', type=float, help='Bin size for discretizing MAF', default=0.01)  # default option does not work
parser.add_argument('compute', metavar='comp', type=int, help='If True, compute quantiles and plots, else runs plotting only', default=1)


argsin = parser.parse_args()
print('\n'.ljust(80, '*'))
print('The following arguments were parsed from file:\n')
print(argsin)
print('\n'.ljust(80, '*'))

pathout = os.path.expanduser(argsin.pathout)
truegt = os.path.expanduser(argsin.truegt)
truegl = os.path.expanduser(argsin.truegl)
pooledgl = os.path.expanduser(argsin.pooledgl)
imputed_gtgp1 = os.path.expanduser(argsin.imp1)  # Beagle
imputed_gtgp2 = os.path.expanduser(argsin.imp2)  # prophaser
datedir = argsin.date
bins_step = argsin.bins
compute = argsin.compute


# compute = True

# Data parameters

x_data = 'maf'
x_bins = np.arange(0.0, 0.5 + bins_step, bins_step) if x_data in ['maf_info', 'maf'] \
            else np.arange(0.0, 1.0 + bins_step, bins_step)
# Custom bins for exact matches counts
x2_bins = [0.0, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.5]  # MAF
lab2_bins = [0.01, 0.03, 0.05, 0.08, 0.15, 0.3, 0.45]  # MAF
lab2_fmt = ['{:.2f}-{:.2f}'.format(i, j) for i, j in zip(x2_bins[:-1], x2_bins[1:])]


# Configure data/plots paths

# pathout = os.path.expanduser('~/PoolImpHuman/results')
# datedir = '20210320'
outdir = os.path.join(pathout, datedir)  # outdir = os.path.join(os.path.expanduser('~/PoolImpHuman/results), datedir)
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Read files

# # pooled HD scenario
# truegt = os.path.expanduser('~/PoolImpHuman/data/20210320/IMP.chr20.snps.gt.vcf.gz')  # '../examples/IMP.chr20.snps.gt.vcf.gz'
# truegl = os.path.expanduser('~/PoolImpHuman/data/20210320/IMP.chr20.snps.gl.vcf.gz')  # '../examples/IMP.chr20.snps.gl.vcf.gz'
# # pooled can also be the file with full LD and missing HD
# # pooledgt = None  # '~/PoolImpHuman/data/20201029/sHG01063.IMP.chr20.missingHD.fullLD.snps.gt.vcf.gz' # Deprecated
# pooledgl = os.path.expanduser('~/PoolImpHuman/data/20210320/IMP.chr20.pooled.snps.gl.vcf.gz')  # '../examples/IMP.chr20.pooled.snps.gl.vcf.gz'
# # imputation with Beagle or with prophaser
# imputed_gtgp1 = os.path.expanduser('~/PoolImpHuman/data/20210320/IMP.chr20.pooled.imputed.vcf.gz')  # '../examples/IMP.chr20.pooled.imputed.vcf.gz'  # prophaser
# imputed_gtgp2 = os.path.expanduser('~/PoolImpHuman/data/20200710-rep/IMP.chr20.pooled.imputed.vcf.gz')  # '../examples/IMP.chr20.pooled.imputed.vcf.gz'  # Beagle


# Build Quality and DataFrame objects for analysis

startT = timeit.default_timer()

dftruegt = vcfdf.PandasMixedVCF(truegt, format='GT', indextype='chrom:pos')

quality1gt = qual.QualityGT(truegt, imputed_gtgp1, 0, idx='chrom:pos')
quality1gl = qual.QualityGL(truegl, imputed_gtgp1, 0, idx='chrom:pos')

if imputed_gtgp2 is not None:
    quality2gt = qual.QualityGT(truegt, imputed_gtgp2, 0, idx='chrom:pos')
    quality2gl = qual.QualityGL(truegl, imputed_gtgp2, 0, idx='chrom:pos')

# Deprecated
# if pooledgt is not None:
#     dfpooled = vcfdf.PandasMixedVCF(pooledgt, format='GT', indextype='chrom:pos')

if pooledgl is not None:
    dfpooled = vcfdf.PandasMixedVCF(pooledgl, format='GL', indextype='chrom:pos')


if compute:

    # Create single-column dataframe for assigning each variant to its MAF-bin

    dX = quality1gt.trueobj.maf
    cutbins = pd.cut(dX.values.squeeze(), x2_bins, labels=lab2_fmt, include_lowest=True)  # lab2_bins
    dBinsX = pd.DataFrame(cutbins, index=dX.index, columns=['binned_' + x_data])
    dBinCounts = dBinsX.reset_index().groupby(by='binned_' + x_data).count()
    print('\r\nNumber of markers per bin:')
    print(dBinCounts)

    totNumMarkers = dftruegt.trinary_encoding().shape[0]
    totNumGenotypes = dftruegt.trinary_encoding().size

    # Filter/mask markers on mismatches and decoding status
    disc1 = quality1gt.diff().dropna()
    miss1 = dfpooled.trinary_encoding()

    disc2 = quality2gt.diff().dropna()  # Non NaN i.e. everything imputed
    miss2 = dfpooled.trinary_encoding()

    unassayed_disc1 = 0.5 * disc1.where(miss1 == -1, other=np.nan)  # keep (mis)matches of only not decoded markers

    unassayed_disc2 = 0.5 * disc2.where(miss2 == -1, other=np.nan)  # discordance is 0, 1, or 2

    # Counts of markers known before imputation, share of the total number of markers over all bins

    decoded = miss1.where(miss1 != -1, other=np.nan)
    dfBinDecoded = dBinsX.join(decoded).groupby(by='binned_' + x_data).count().mean(axis=1) # / totNumMarkers
    dfBinDecoded.name = 'fully_decoded'
    print('\r\nAverage number of decoded markers per bin over all samples:')
    print(dfBinDecoded)

    # Additional exact matches after imputation

    # Count not decoded and correctly imputed
    unassayed_impok1 = unassayed_disc1.where(unassayed_disc1 == 0.0, other=np.nan)  # keep not decoded and full matches after imputation
    counts_unassayed_impok1 = unassayed_impok1.notnull().mean(axis=1)  # mean overall sample per variant
    counts_unassayed_impok1.name = 'counts'

    unassayed_impok2 = unassayed_disc2.where(unassayed_disc2 == 0.0, other=np.nan)  # keep not decoded and full matches after imputation
    counts_unassayed_impok2 = unassayed_impok2.notnull().mean(axis=1)
    counts_unassayed_impok2.name = 'counts'

    # Count not decoded per variant
    not_decoded1 = miss1.where(miss1 == -1, other=np.nan)
    counts_not_decoded1 = not_decoded1.notnull().mean(axis=1)  # mean overall sample per variant

    not_decoded2 = miss2.where(miss2 == -1, other=np.nan)
    counts_not_decoded2 = not_decoded2.notnull().mean(axis=1)

    # Proportion of full matches after imputation for not decoded markers

    sExact1 = counts_unassayed_impok1 # / counts_not_decoded1
    sExact1.name = 'exact_matches'  # counts exactly imputed only
    # exact_matches1 = qual.QuantilesDataFrame(dX, sExact1)
    dfmatch1 = dBinsX.join(sExact1).groupby(by='binned_' + x_data).sum()
    dfmatch1['dataset'] = '1'

    sExact2 = counts_unassayed_impok2 #/ counts_not_decoded2
    sExact2.name = 'exact_matches'  # counts exactly imputed only
    # exact_matches2 = qual.QuantilesDataFrame(dX, sExact2)
    dfmatch2 = dBinsX.join(sExact2).groupby(by='binned_' + x_data).sum()
    dfmatch2['dataset'] = '2'

    dfmatches = pd.concat([dfmatch1, dfmatch2])
    print('\r\nAverage number of markers that are additionally imputed correctly over all samples:')
    print(dfmatches)

    # Total of exact matches after pooling and imputation, proportion per bin and out of total markers

    dfTotMatch = dfmatches.join(dfBinDecoded)
    dfTotMatch['exact_matches'] = dfTotMatch['exact_matches'] + dfTotMatch['fully_decoded']
    dfTotMatch['prop_exact_geno'] = dfTotMatch['exact_matches'] / dBinCounts['variants']
    dfTotMatch['prop_fully_deco'] = dfTotMatch['fully_decoded'] / dBinCounts['variants']
    dfTotMatch = dfTotMatch[['dataset', 'fully_decoded', 'exact_matches', 'prop_fully_deco', 'prop_exact_geno']]  # reorder columns

    tabmatches = dfTotMatch.reset_index().pivot(index='binned_' + x_data, columns='dataset').T
    print('\r\nAverage number of markers that are correctly genotyped over all samples:')
    print(tabmatches)
    tabmatches.to_latex(buf=os.path.join(outdir, 'tables-decoded-imputed-bin.tex'),
                        sparsify=True,
                        multirow=True,
                        float_format="%.3f",
                        caption='Number of genotyped markers after pooling and after imputation per data {} bin'.format(
                            'MAF' if (x_data == 'maf' or x_data == 'maf_info') else 'AAF'),
                        label='tab:tables-decoded-imputed-bin')

