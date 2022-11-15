import sys, os
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
print(rootdir)
sys.path.insert(0, rootdir)
import pandas as pd
import numpy as np

from genotypooler.poolSNPs.metrics import quality
import subprocess
import argparse
import matplotlib.pyplot as plt


"""
Compute results with customized metrics from true vs. imputed data sets

Usage:
$ python3 -u imputation_quality.py <path to directory> <VCF file with true genotypes> <VCF file with imputed genotypes> <path to script converting GT to GL> <ID type for variants>

Ex.
$ python3 -u  imputation_quality.py ./ IMP.chr20.snps.gt.vcf.gz IMP.chr20.pooled.imputed.vcf.gz ../bin/gt_to_gl.sh 'id'
"""

parser = argparse.ArgumentParser(description='Compute and plot'
                                             'customized imputation accuracy metrics')
parser.add_argument('directory', metavar='dir', type=str, help='Path to directory with files', default=None)
parser.add_argument('true', metavar='tru', type=str, help='File with true genotypes', default=None)
parser.add_argument('imputed', metavar='imp', type=str, help='Imputed file with genotypes (GT:DS:GP)', default=None)
parser.add_argument('gconverter', metavar='gcv', type=str, help='Path to script converting GT to GL', default='~/PoolImpHuman/bin/bash-src/gt_to_gl.sh')
parser.add_argument('ident', metavar='idt', type=str, help='Type of identifier for variants (id/chrom:pos)', default='id')

argsin = parser.parse_args()
dirin = argsin.directory
ftrue = argsin.true
fimp = argsin.imputed
gconv = argsin.gconverter
idt = argsin.ident

paths = {'beaglegt': {
    'true': os.path.join(dirin, ftrue),
    'imputed': os.path.join(dirin, fimp)},
    'beaglegl': {
        'true': os.path.join(dirin, ftrue.replace('.gt.', '.gl.')),
        'imputed': os.path.join(dirin, fimp)},
}

convertgtgl = True
if convertgtgl:
    cmd = 'bash {} {} {}'.format(gconv, paths['beaglegt']['true'], paths['beaglegl']['true'])
    subprocess.run(cmd, shell=True,)

qbeaglegt = quality.QualityGT(*paths['beaglegt'].values(), 0, idx=idt)
qbeaglegl = quality.QualityGL(paths['beaglegl']['true'], paths['beaglegl']['imputed'], 0, idx=idt)

try:
    entro = qbeaglegl.cross_entropy
except KeyError:
    entro = None

tabbeaglegl = pd.concat([qbeaglegt.concordance(),
                         qbeaglegt.trueobj.af_info,
                         qbeaglegt.pearsoncorrelation(),
                         qbeaglegt.precision,
                         qbeaglegt.accuracy,
                         qbeaglegt.recall,
                         qbeaglegt.f1_score,
                         entro], axis=1)
dosbeaglegl = qbeaglegt.alleledosage()

tabbeaglegl.head()

plt.rcParams["figure.figsize"] = [5*4, 4*2]
fig, axes = plt.subplots(2, 4)

tabbeaglegl.plot.scatter('af_info', 'precision_score', ax=axes[0, 0], s=0.7)
axes[0, 0].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'accuracy_score', ax=axes[0, 1], s=0.7)
axes[0, 1].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'concordance', ax=axes[0, 2], s=0.7)
axes[0, 2].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'f1_score', ax=axes[0, 3], s=0.7)
axes[0, 3].set_ylim(0.0, 1.0)
tabbeaglegl.plot.scatter('af_info', 'r_squared', ax=axes[1, 0], s=0.7)
axes[1, 0].set_ylim(-0.2, 1.0)
if entro is not None:
    tabbeaglegl.plot.scatter('af_info', 'cross_entropy', ax=axes[1, 1], s=0.7)
    axes[1, 1].set_ylim(-0.5, 5.0)
axes[1, 2].scatter(dosbeaglegl[0], dosbeaglegl[1], s=0.7)
axes[1, 2].set_xlabel('true allele dosage')
axes[1, 2].set_ylabel('imputed allele dosage')
axes[1, 2].set_ylim(0.0, 2.0)

for ax in axes.flatten()[:-2]:
    # cast color to white 'w' if dark background
    ax.set_xlabel('true alternate allele frequency', color='k')
    ax.set_ylabel(ax.get_ylabel().replace('_', ' '), color='k')
plt.suptitle('Evaluation of imputation performance')
plt.savefig(os.path.join(os.path.dirname(paths['beaglegt']['imputed']), 'imputation_quality_gtgl.pdf'))
plt.show()
