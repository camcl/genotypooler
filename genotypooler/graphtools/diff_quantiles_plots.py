"""
Plots comparing the accuracy of imputation between consecutive runs.
The following metrics are computed across MAF range and shown with quantiles [25%, 50%, 75%] dispersion.
* Concordance: marker-wise difference cycleN - cycleN+1
* Cross-entropy: marker-wise difference cycleN - cycleN+1

Quantiles computed on rolling windows.

Arguments are parsed from a file where they must be written in this order:
pathout <>
truegt <>
truegl <>
imp1 <path/to/file/imputed/with/data1>
imp2 <path/to/file/imputed/with/data2>
date <>
rollwin <>
bins <>
compute <>

Beware that strings for paths should be written just as text (without quotes!) in the argsfile!

Command line usage (assuming the current directory is genotypooler/graphtools)
$ python3 -u diff_quantiles_plots.py @/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/cycle2/quantilesargsfile.txt
"""

import os, sys
import numpy as np
import pandas as pd
import json
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

import warnings
warnings.filterwarnings("ignore")

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs.metrics import quality as qual
from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.persotools.files import *


### COMMAND-LINE PARSING AND PARAMETERS (arguments are parsed from a file
parser = argparse.ArgumentParser(description='Plots the imputation performance'
                                             '(concordance and cross-entropy)',
                                 fromfile_prefix_chars='@')
parser.add_argument('pathout', metavar='out', type=str, help='Results directory', default=None)
parser.add_argument('truegt', metavar='gt', type=str, help='File with true data (GT format)', default=None)
parser.add_argument('truegl', metavar='gl', type=str, help='File with true data (GL format)', default=None)
parser.add_argument('imp1', metavar='imp1', type=str, help='File with imputed data 1 (GT and/or GP formats)', default=None)
parser.add_argument('imp2', metavar='imp2', type=str, help='File with imputed data 2 (GT and/or GP formats)', default=None)
parser.add_argument('date', metavar='date', type=str, help='Date of the experiment (YYYMMDD)', default=None)
parser.add_argument('rollwin', metavar='wq', type=int, help='Number of markers per rolling window', default=1000)  # default option does not work
parser.add_argument('bins', metavar='bin', type=float, help='Bin size for discretizing MAF', default=0.01)  # default option does not work
parser.add_argument('compute', metavar='comp', type=int, help='If True, compute quantiles and plots, else runs plotting only', default=1)
parser.add_argument('mask', metavar='mask', type=str, help='File with filtered data (GP format)', default=None)
parser.add_argument('changed', metavar='chg', type=int, help='If True, mask unchanged priors', default=None)

argsin = parser.parse_args()
print('\n'.ljust(80, '*'))
print('The following arguments were parsed from file:\n')
print(argsin)
print('\n'.ljust(80, '*'))

pathout = argsin.pathout
truegt = argsin.truegt
truegl = argsin.truegl
imputed_1 = argsin.imp1
imputed_2 = argsin.imp2
datedir = argsin.date
rQ = argsin.rollwin * 4  # try some integers in the range [2, 5] and see what is the smoothest
bS = argsin.bins
compute = argsin.compute
fmask = argsin.mask
changed = argsin.changed

# Create gbool mask for data

# parent_dir = "/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/"
# vcf_changes = parent_dir + "cycle6/STU.Chr1.SNPs.pruned.sorted.pooled.changed_priors.vcf"
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

# Data parameters

x_data = 'binned_maf'

# Configure data/plots paths

outdir = os.path.join(pathout, datedir)
if not os.path.exists(outdir):
    os.mkdir(outdir)


print('\r\nData written to {}'.format(outdir))

# Plot styling

# Specific to this plotting script
sns.set(rc={'figure.figsize': (10, 8)})
sns.set_style('whitegrid')
titlesz = 24
axlabsz= 20
axticksz = 16
legsz = 20
yscale = {
    'concordance': (-0.05, 0.01),
    'cross_entropy': (-0.15, 0.45)
}


# Function/Tools

def diff_rollquants(dX: pd.DataFrame, dS1: pd.Series) -> pd.DataFrame:
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bS)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = ['cycleN - cycleN+1'] * pctY1.shape[0]
    diff_rollquants = pctY1

    return diff_rollquants


# Load data and check

q1gt = qual.QualityGT(truegt, imputed_1, 0, idx='chrom:pos', mask=gbool)
q1gl = qual.QualityGL(truegl, imputed_1, 0, idx='chrom:pos', mask=gbool)
if fmask != '.' and changed:  # parse the complementary data i.e. mask changed priors
    q12gt = qual.QualityGT(truegt, imputed_1, 0, idx='chrom:pos', mask=~gbool)
    q12gl = qual.QualityGL(truegl, imputed_1, 0, idx='chrom:pos', mask=~gbool)

print('\r\n{} variants from {} samples read from {}'.format(len(q1gt.trueobj.variants),
                                                            len(q1gt.trueobj.samples),
                                                            os.path.basename(truegt)))
print('\r\n{} variants from {} samples read from {}'.format(len(q1gt.imputedobj.variants),
                                                            len(q1gt.imputedobj.samples),
                                                            os.path.basename(imputed_1)))
if compute:
    bgldiff = q1gt.diff()

q2gt = qual.QualityGT(truegt, imputed_2, 0, idx='chrom:pos', mask=gbool)
q2gl = qual.QualityGL(truegl, imputed_2, 0, idx='chrom:pos', mask=gbool)
if fmask != '.' and changed:  # parse the complementary data i.e. mask changed priors
    q22gt = qual.QualityGT(truegt, imputed_2, 0, idx='chrom:pos', mask=~gbool)
    q22gl = qual.QualityGL(truegl, imputed_2, 0, idx='chrom:pos', mask=~gbool)

print('\r\n{} variants from {} samples read from {}'.format(len(q2gt.trueobj.variants),
                                                            len(q2gt.trueobj.samples),
                                                            os.path.basename(truegl)))
print('\r\n{} variants from {} samples read from {}'.format(len(q2gl.imputedobj.variants),
                                                            len(q2gl.imputedobj.samples),
                                                            os.path.basename(imputed_2)))

mafS = q1gt.trueobj.maf

if compute:
    metrics = {
        'concordance': {'1': q1gt.marker_wise_concordance(),
                        '2': q2gt.marker_wise_concordance()},
        'cross_entropy': {'1': q1gl.marker_wise_cross_entropy,
                          '2': q2gl.marker_wise_cross_entropy}
    }
    if fmask != '.' and changed:  # parse the complementary data i.e. mask changed priors
        antimetrics = {
            'concordance': {'1': q12gt.marker_wise_concordance(),
                            '2': q22gt.marker_wise_concordance()},
            'cross_entropy': {'1': q12gl.marker_wise_cross_entropy,
                              '2': q22gl.marker_wise_cross_entropy}
        }

dataquants = {
    'concordance': os.path.join(outdir, 'diff_rolling_quantiles_concordance.json'),
    'cross_entropy': os.path.join(outdir, 'diff_rolling_quantiles_cross_entropy.json')
}
if fmask != '.' and changed:
    antidataq = {
        'concordance': None,
        'cross_entropy': None
    }

# Process and write data

if compute:
    for metric, d in metrics.items():
        if d is not None:
            yDF_1, yDF_2 = list(d.values())
            yS1_2 = yDF_1.subtract(yDF_2).mean(axis=1).rename(metric)

            # Compute quantiles
            print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
            pctY_comp = diff_rollquants(mafS, yS1_2)

            # Compute mean over all markers
            print('Computing means for {}'.format(metric).ljust(80, '.'))
            pctY_comp['mean'] = pctY_comp['dataset'].apply(lambda x: yS1_2.mean())
            jsonf = dataquants[metric]
            pctY_comp.to_json(jsonf,
                              orient='records')

    if fmask != '.' and changed:  # parse the complementary data i.e. mask changed priors
        for metric, d in antimetrics.items():
            if d is not None:
                zDF_1, zDF_2 = list(d.values())
                zS1_2 = zDF_1.subtract(zDF_2).mean(axis=1).rename(metric)
                # Compute quantiles
                print('Computing quantiles for {} in complementary data'.format(metric).ljust(80, '.'))
                pctZ_comp = diff_rollquants(mafS, zS1_2)
                # Compute mean over all markers
                print('Computing means for {} in complementary data'.format(metric).ljust(80, '.'))
                pctZ_comp['mean'] = pctZ_comp['dataset'].apply(lambda x: zS1_2.mean())
                antidataq[metric] = pctZ_comp

# Read processed reshaped data for plotting and draw figures

sns.set(font_scale=1.75)  # multiplication factor!
sns_palette = "husl"  # Choose a palette with  fixed number of levels, or fix n_colors with husl

for dquant, f in dataquants.items():
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        meanf = {}

        if fmask != '.' and changed:
            dataf['dataset'] = 'changed priors'
            antidataf = antidataq[dquant]
            antidataf['dataset'] = 'unchanged priors'
            dataf = pd.concat([dataf, antidataf], axis=0).reset_index(drop=True)
        print(dataf)

        gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x=x_data, y=dquant,
                          hue='dataset', palette=sns.color_palette(sns_palette, n_colors=2), linewidth=1)
        for i, dset in enumerate(np.unique(dataf['dataset'].values)):  # enumerate(['cycleN - cycleN+1']):
            df = dataf[dataf['dataset'] == dset]
            meanf[dset] = df['mean'].mean()
            gY.fill_between(df[df.quantiles == 1.0][x_data],
                            df[df.quantiles == 0.0][dquant],
                            df[df.quantiles == 1.0][dquant],
                            color=sns.color_palette(sns_palette, n_colors=2)[i],
                            alpha=0.1)
            gY.fill_between(df[df.quantiles == 0.99][x_data],
                            df[df.quantiles == 0.01][dquant],
                            df[df.quantiles == 0.99][dquant],
                            color=sns.color_palette(sns_palette, n_colors=2)[i],
                            alpha=0.25)
            gY.fill_between(df[df.quantiles == 0.75][x_data],
                            df[df.quantiles == 0.25][dquant],
                            df[df.quantiles == 0.75][dquant],
                            color=sns.color_palette(sns_palette, n_colors=2)[i],
                            alpha=0.40)
        gY.set_xlabel('True minor allele frequency in {} population'.format('study' if x_data == 'binned_maf'
                                                                            else 'main'),
                      fontsize=axlabsz)
        gY.set_ylabel(f"Difference in {dataf.columns[2].replace('_', ' ')} (cycleN - cycleN+1)", fontsize=axlabsz)
        gY.set(ylim=yscale[dquant])
        if gbool is not None and not (fmask != '.' and changed):
            gY.set_title(f'''Number of genotypes used (variants x samples) = {gbool.sum()} 
            i.e. {gbool.sum() * 100 / gbool.size:2.1f}% of the dataset''', fontsize=16)

        handles, labels = gY.get_legend_handles_labels()
        if fmask != '.' and changed:
            labels[-2] = '{} (mean = {:.5f})'.format(labels[-2], meanf[np.unique(dataf['dataset'].values)[0]])
        labels[-1] = '{} (mean = {:.5f})'.format(labels[-1], meanf[np.unique(dataf['dataset'].values)[-1]])  # meanf['cycleN - cycleN+1'])
        gY.legend(handles, labels, loc='lower right' if dquant == 'concordance' else 'upper right', fontsize=legsz)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'diff_{}_percentiles_rQ={}_bS={}_xdata={}.pdf'.format(dquant, rQ, bS, x_data.lstrip('binned_'))))
        # plt.show()
        plt.close()