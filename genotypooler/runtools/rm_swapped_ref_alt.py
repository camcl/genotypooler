import os, sys
import argparse
import subprocess
import numpy as np

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import dataframe as vcfdf
from genotypooler.poolSNPs import pybcf
from genotypooler.persotools.files import *

'''
Removes the markers where the reference and the alternate alleles are opposite in the VCF files
for the founders and for the inbred lines. 
'''

CHROM_TABLE = {'1A': '1',
               '1B': '2',
               '1D': '3',
               '2A': '4',
               '2B': '5',
               '2D': '6',
               '3A': '7',
               '3B': '8',
               '3D': '9',
               '4A': '10',
               '4B': '11',
               '4D': '12',
               '5A': '13',
               '5B': '14',
               '5D': '15',
               '6A': '16',
               '6B': '17',
               '6D': '18',
               '7A': '19',
               '7B': '20',
               '7D': '21'}

parser = argparse.ArgumentParser(description='Remove loci with swapped alleles.')
parser.add_argument('indir1', metavar='indir1', type=str, help='Directory with the data 1 to process', default=None)
parser.add_argument('indir2', metavar='indir2', type=str, help='Directory with the data 2 to process', default=None)
parser.add_argument('founders', metavar='founders', type=str, help='VCF file with data 1 (GT format)', default=None)
parser.add_argument('inbred', metavar='inbred', type=str, help='VCF file with data 1 (GT format)', default=None)


argsin = parser.parse_args()
print('\n'.ljust(80, '*'))
print('The following arguments were parsed from the command line:\n')
print(argsin)
print('\n'.ljust(80, '*'))

dir_in1 = argsin.indir1  # '/home/camille/MagicWheat/src/genotypooler/examples'
dir_in2 = argsin.indir2
foundersfile = argsin.founders  # 'Chr1.Founders.vcf.gz'
inbredfile = argsin.inbred  # 'Chr1.SNPs.pruned.nomiss.vcf.gz'
dir_all = 'unfiltered-opposite-markers'


def find_swap(x):
    """Works on 1D numpy slices: if (x['ref_founders'] == x['alt_inbred']
                            and x['ref_inbred'] == x['alt_founders']) else 0"""
    return 1 if (x[0] != x[2]) else np.nan

# Copy original files to a separate folder

for dir_in in [dir_in1, dir_in2]:
    if not os.path.exists(os.path.join(dir_in, dir_all)):
        os.makedirs(os.path.join(dir_in, dir_all))

for dir_in, f_gz in [(dir_in1, foundersfile), (dir_in2, inbredfile)]:
    # break  # comment if original files contain swapped alleles
    cmd_copy = ' '.join(['bcftools',
                         'view -Oz -o',
                         os.path.join(dir_in, dir_all, f_gz),
                         '--threads {}'.format(os.cpu_count()),
                         os.path.join(dir_in, f_gz)
                         ])
    print(cmd_copy)
    subprocess.run(cmd_copy, shell=True)
    pybcf.index(f_gz, os.path.join(dir_in, dir_all))

dffounders = vcfdf.PandasMixedVCF(os.path.join(dir_in1, dir_all, foundersfile), format='GT', indextype='id')

dfinbred = vcfdf.PandasMixedVCF(os.path.join(dir_in2, dir_all, inbredfile), format='GT', indextype='id')

dfalleles = dffounders.alleles.join(dfinbred.alleles, lsuffix='_founders', rsuffix='_inbred')
print(dfalleles)

dfswapped = dfalleles.apply(find_swap, axis=1).dropna()
print(f'\nFound {len(dfswapped.index)} allele swaps at the following markers'.ljust(80, '.'))
print(dfswapped.index)

print('\nCreating VCF without swapped markers'.ljust(80, '.'))
swapped_chrom_pos = dfswapped.index.str.split(':').to_list()
for dir_in, f_gz in [(dir_in1, foundersfile), (dir_in2, inbredfile)]:
    with open(os.path.join(dir_in, dir_all, 'swapped_chrom_pos.coords'), 'w') as fcoord:
        fcoord.writelines(CHROM_TABLE[chrom] + '\t' + pos + '\n' for chrom, pos in swapped_chrom_pos)
    if len(swapped_chrom_pos) > 0:  # found swaps
        cmd_filter = ' '.join(['bcftools',
                               'view -Oz -o',
                               os.path.join(dir_in, f_gz),
                               '-T ^' + os.path.join(dir_in, dir_all, 'swapped_chrom_pos.coords'),  # use "-R ^" if the VCF files are indexed
                               # '--threads {}'.format(os.cpu_count()),
                               os.path.join(dir_in, dir_all, f_gz)
                               ])
        subprocess.run(cmd_filter, shell=True, cwd=os.path.join(dir_in, dir_all))
        pybcf.index(f_gz, dir_in)

# TODO: (next step in the workflow) convert study population from GT to GL
# (venv3.6) (base) camille@camille-Precision-3551:~/MagicWheat/src/genotypooler/examples$ bash ../bin/gt_to_log_gl.sh Chr1.SNPs.pruned.nomiss.vcf.gz Chr1.SNPs.pruned.nomiss.gl.vcf.gz