"""
Plots quantifying the imputation performance across multiple cycles of decoding and imputation.
The following metrics are computed across MAF range and shown with quantiles [25%, 50%, 75%] dispersion.
* Concordance:
* Cross-entropy:

Quantiles computed on rolling windows.

Arguments are parsed from a file where they must be written in this order:
pathout <>
truegt <>
truegl <>
ncycles <list of (integer) cycles to plot, separated by commas without space>
impdir <path/pattern/to/cycle/dir/with/file/imputed/>
impfile <name_file_imputed>
date <>
rollwin <>
bins <>
compute <>
mask <>

Beware that strings for paths should be written just as text (without quotes!) in the argsfile!

Command line usage (assuming the current directory is genotypooler/graphtools)
$ python3 -u multi_quantiles_plots.py @/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/cycle2/quantilesargsfile.txt
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
parser.add_argument('ncycles', metavar='ncyc', type=str, help='File with imputed data 1 (GT and/or GP formats)', default=None)
parser.add_argument('impdir', metavar='impd', type=str, help='File with imputed data 2 (GT and/or GP formats)', default=None)
parser.add_argument('impfile', metavar='impf', type=str, help='File with imputed data 2 (GT and/or GP formats)', default=None)
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
ncycles = argsin.ncycles.split(',')
impdir = argsin.impdir
impfile = argsin.impfile
datedir = argsin.date
rQ = argsin.rollwin
bS = argsin.bins
compute = argsin.compute
fmask = argsin.mask
changed = argsin.changed

# No gbool mask for data

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
sns.set(rc={'figure.figsize': (10, 8)})  # specific to this plotting sripts
sns.set_style('whitegrid')
titlesz = 24
axlabsz= 20
axticksz = 16
legsz = 20
yscale = {
    'concordance': (0.75, 1.0),  # no zoom-in: (0.0, 1.0),
    'cross_entropy': (0.0, 1.5)  # no zoom-in: (0.0, 12.0)
}


# Function/Tools

def single_rollquants(dX: pd.DataFrame, dS1: pd.Series, N: int) -> pd.DataFrame:
    pdf1 = qual.QuantilesDataFrame(dX,
                                   dS1,
                                   bins_step=bS)
    pctY1 = pdf1.binnedX_rolling_quantilY(rollwin=rQ)
    pctY1['dataset'] = [str(N)] * pctY1.shape[0]
    single_rollquants = pctY1

    return single_rollquants


# Load data for first (default) cycle

q1gt = qual.QualityGT(truegt, os.path.join(impdir + ncycles[0], impfile), 0, idx='chrom:pos', mask=gbool)
q1gl = qual.QualityGL(truegl, os.path.join(impdir + ncycles[0], impfile), 0, idx='chrom:pos', mask=gbool)

print('\r\n{} variants from {} samples read from {}'.format(len(q1gt.trueobj.variants),
                                                            len(q1gt.trueobj.samples),
                                                            os.path.basename(truegt)))
print('\r\n{} variants from {} samples read from {}'.format(len(q1gt.imputedobj.variants),
                                                            len(q1gt.imputedobj.samples),
                                                            impdir + ncycles[0]))
if compute:
    bgldiff = q1gt.diff()
    metrics = {
        'concordance': {'1': q1gt.concordance()},
        'cross_entropy': {'1': q1gl.cross_entropy}
    }

mafS = q1gt.trueobj.maf

#TODO: Loop over other imputed datasets in the other cycles, append to metrics with key=cycle num
for n in ncycles[1:]:
    q2gt = qual.QualityGT(truegt, os.path.join(impdir + n, impfile), 0, idx='chrom:pos', mask=gbool)
    q2gl = qual.QualityGL(truegl, os.path.join(impdir + n, impfile), 0, idx='chrom:pos', mask=gbool)

    print('\r\n{} variants from {} samples read from {}'.format(len(q2gt.trueobj.variants),
                                                                len(q2gt.trueobj.samples),
                                                                os.path.basename(truegl)))
    print('\r\n{} variants from {} samples read from {}'.format(len(q2gl.imputedobj.variants),
                                                                len(q2gl.imputedobj.samples),
                                                                impdir + n))
    if compute:
        metrics['concordance'][n] = q2gt.concordance()
        metrics['cross_entropy'][n] = q2gl.cross_entropy


dataquants = {
    'concordance': os.path.join(outdir, 'multi_rolling_quantiles_concordance.json'),
    'cross_entropy': os.path.join(outdir, 'multi_rolling_quantiles_cross_entropy.json')
}

# Process and write data

if compute:
    for metric, d in metrics.items():
        print('Computing quantiles for {}'.format(metric).ljust(80, '.'))
        # initialize with data from cycle 1
        yS_1 = d['1']
        pctY_multi = single_rollquants(mafS, yS_1, 1)
        print('Computing means for {}'.format(metric).ljust(80, '.'))
        pctY_multi['mean'] = pctY_multi['dataset'].apply(lambda x: yS_1.mean())

        for n in ncycles[1:]:
            yS_N = d[n]
            # Compute quantiles
            pctY_N = single_rollquants(mafS, yS_N, int(n))
            # Compute mean over all markers
            pctY_N['mean'] = pctY_N['dataset'].apply(lambda x: yS_N.mean())
            pctY_multi = pd.concat([pctY_multi, pctY_N], axis=0).reset_index(drop=True)

        jsonf = dataquants[metric]
        pctY_multi.to_json(jsonf, orient='records')

# Read processed reshaped data for plotting and draw figures

sns.set(font_scale=1.75)  # multiplication factor!

for dquant, f in dataquants.items():
    if f is not None:
        dataf = pd.read_json(f, orient='records')
        print(dataf)
        meanf = {}

        gY = sns.lineplot(data=dataf[dataf.quantiles == 0.5], x=x_data, y=dquant,
                          hue='dataset', palette="colorblind", linewidth=1)
        for i, dset in enumerate(ncycles):
            df = dataf[dataf['dataset'] == int(dset)]
            meanf[dset] = df['mean'].mean()
            if i == 0 or i == len(ncycles) - 1:
                gY.fill_between(df[df.quantiles == 1.0][x_data],
                                df[df.quantiles == 0.0][dquant],
                                df[df.quantiles == 1.0][dquant],
                                color=sns.color_palette('colorblind')[i],
                                alpha=0.05)
                gY.fill_between(df[df.quantiles == 0.99][x_data],
                                df[df.quantiles == 0.01][dquant],
                                df[df.quantiles == 0.99][dquant],
                                color=sns.color_palette('colorblind')[i],
                                alpha=0.15)
                gY.fill_between(df[df.quantiles == 0.75][x_data],
                                df[df.quantiles == 0.25][dquant],
                                df[df.quantiles == 0.75][dquant],
                                color=sns.color_palette('colorblind')[i],
                                alpha=0.25)
        gY.set_xlabel('True minor allele frequency in {} population'.format('study' if x_data == 'binned_maf'
                                                                            else 'main'),
                      fontsize=axlabsz)
        gY.set_ylabel(str.capitalize(dataf.columns[2].replace('_', ' ')), fontsize=axlabsz)
        gY.set(ylim=yscale[dquant])
        handles, labels = gY.get_legend_handles_labels()
        for i, leglab in enumerate(labels):
            labels[i] = 'Cycle {} (mean = {:.5f})'.format(leglab, meanf[ncycles[i]])
        gY.legend(handles, labels, loc='lower left' if dquant == 'concordance' else 'upper right', fontsize=legsz)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, 'multi_{}_percentiles_rQ={}_bS={}_xdata={}.pdf'.format(dquant, rQ, bS, x_data.lstrip('binned_'))))
        plt.close()
