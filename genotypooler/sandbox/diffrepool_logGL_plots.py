"""
Plots comparing ...
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

last_cycle = 42
outdir = "/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/20230608"
otherdirs = []
# ["/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/20230605",
#              "/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/20230608"]
wfactor = 0.01
pooledf = "STU.Chr1.SNPs.pruned.sorted.pooled.vcf.gz"
jsonout = "diff_pooled_log_gl.json"
figout = "diff_pooled_log_gl.pdf"
compute = 0

fmask = '.'  # '/home/camille/IterDecodeImpute/runs/poolimputeSNPs/results/data/1/20230605/cycle2/STU.Chr1.SNPs.pruned.sorted.pooled.changed_priors.vcf.gz'  # argsin.mask
changed = 1  # argsin.changed

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

# Plot styling

# Specific to this plotting script
sns.set(rc={'figure.figsize': (12, 8)})
sns.set_style('whitegrid')
sns.set(font_scale=1.75)  # multiplication factor!
titlesz = 24
axlabsz= 20
axticksz = 16
legsz = 20


# Function/Tools

def log_geno_tuples_diff(sGtup1: pd.Series, sGtup2: pd.Series) -> pd.Series:
    dG1 = pd.DataFrame.from_records(sGtup1.values,
                                    index=sGtup1.index,
                                    columns=['RR', 'RA', 'AA']).astype(float)
    dG2 = pd.DataFrame.from_records(sGtup2.values,
                                    index=sGtup2.index,
                                    columns=['RR', 'RA', 'AA']).astype(float)

    return dG2.sub(dG1).sum(axis=1)


# Load data and check

for n in range(2, last_cycle + 1):
    if not compute:
        break
    else:
        loggl1 = vcfdf.PandasMixedVCF(os.path.join(outdir, f'cycle{n - 1}', pooledf),
                                      format='GL',
                                      indextype='chrom:pos')
        loggl2 = vcfdf.PandasMixedVCF(os.path.join(outdir, f'cycle{n}', pooledf),
                                      format='GL',
                                      indextype='chrom:pos')

        loggldiff = loggl1.genotypes().combine(loggl2.genotypes(), log_geno_tuples_diff)
        loggldiff.to_json(os.path.join(outdir, f'cycle{n}', jsonout), orient='records')

# Printouts for verification
if compute:
    print(loggl1.genotypes())
    print('\r\n')
    # find a way to subtract GL tuple by GL tuple
    print(log_geno_tuples_diff(loggl1.genotypes()['A1310_A1310'], loggl2.genotypes()['A1310_A1310']))
    print('\r\n')

per_cycle_diff = pd.DataFrame(data=np.nan,
                              index=range(1, last_cycle),
                              columns=['tot_log_gl_diff'])  # not relevant for the initial cycle
per_cycle_diff.index.name = 'cycle'
for n in range(2, last_cycle + 1):
    loggldiff = pd.read_json(os.path.join(outdir, f'cycle{n}', jsonout), orient='records')
    per_cycle_diff.loc[n] = loggldiff.sum().sum()
per_cycle_diff.reset_index(drop=False, inplace=True)
per_cycle_diff['w'] = wfactor
per_cycle_diff.to_json(os.path.join(outdir, jsonout), orient='records')

for otherf in otherdirs:
    try:
        dfotherw = pd.read_json(os.path.join(otherf, jsonout), orient='records')
        per_cycle_diff = pd.concat([per_cycle_diff, dfotherw], axis=0, ignore_index=True)
    except (ValueError, FileNotFoundError):
        print('No data available')

print('\r\n')
print(per_cycle_diff)

gLogGLDiff = sns.lineplot(data=per_cycle_diff, x='cycle', y='tot_log_gl_diff',
                          hue='w', palette="viridis", linewidth=3, legend="auto")
gLogGLDiff.set_xlabel("Cycle N")
gLogGLDiff.set_ylabel("Total difference in log-GL")
gLogGLDiff.set_title("Changes in the genotype priors in the decoded data between cycles N and N-1")
for outd in [outdir, *otherdirs]:
    plt.savefig(os.path.join(outd, figout))

plt.yscale('log')
for outd in [outdir, *otherdirs]:
    try:
        plt.savefig(os.path.join(outd, 'log_yscale_' + figout))
    except:
        print('No data directory')
