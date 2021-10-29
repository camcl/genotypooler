import sys, os
import numpy as np
import pandas as pd
from itertools import starmap, repeat
import shutil
import subprocess
import multiprocessing as mp
import argparse
import timeit
import datetime

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, './')

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
print(rootdir)
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import poolvcf
from genotypooler.poolSNPs import pybcf

from persotools.files import delete_file, mkdir

'''
Parallelized file processing
For VCF-file bigger than some dozen of thousands of varaiants (e.g. 1 million SNPs).
Can be run on Rackham at Uppmax.

Steps:
* Read main vcf and write chunks: bash script based on bcftools
* Simulate pooling on chunks and decode into GT
* Merge pooled chunks back to read-for-use compressed VCF file

Usage: 
$ python3 parallel_pooling.py data/IMP.chr20.snps.gt.vcf.gz data/IMP.chr20.pooled.snps.gt.vcf.gz 4
'''

### COMMAND-LINE PARSING AND PARAMETERS
parser = argparse.ArgumentParser(description='Run parallelized pooling simulation'
                                             'on the whole set of samples')
parser.add_argument('pathin', metavar='in', type=str, help='File to pool', default=None)
parser.add_argument('pathout', metavar='out', type=str, help='Pooled file', default=None)
parser.add_argument('cores', metavar='nc', type=int, nargs='?', help='Number of cores to use', default=None)

argsin = parser.parse_args()

nb_cores = os.cpu_count() if argsin.cores is None else argsin.cores
try:
    tmp_path = os.environ['SNIC_TMP']
except KeyError:
    tmp_path = os.path.join(proj_dir, 'data/tmp')

try:
    assert os.path.exists('snic')
    snic_proj = '/crex/proj/snic2019-8-216/private'
except AssertionError:
    snic_proj = os.path.join(proj_dir, 'data')
print('SNIC PROJ: {}'.format(snic_proj))

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Input file = {}'.format(os.path.expanduser(argsin.pathin)))
print('Output file = {}'.format(os.path.expanduser(argsin.pathout)))
print('\n'.rjust(80, '*'))

### SIMULATE POOLING ON PACKED DATA CHUNKS
data_dir = './data'
os.chdir(tmp_path)
fingz = os.path.expanduser(argsin.pathin)  # /home/camille/...
basin = os.path.basename(fingz).rstrip('.gz')
basout = os.path.basename(argsin.pathout)


start = timeit.default_timer()
prename = 'pack'

wd = os.path.join(snic_proj, 'tmp')
files0 = [f for f in os.listdir(wd)]
# keep only chunks to process
for f0 in files0:
    if not f0.startswith(prename):
        files0.remove(f0)

indices = np.arange(len(files0))

print('\r\n{} files found will be pooled'.format(len(files0)).ljust(80, '.'))
print([os.path.join(wd, f0) for f0 in files0])

files1 = ['pooled.{}'.format(f0).rstrip('.gz') for f0 in files0]
args1 = list(zip([os.path.join(wd, f0) for f0 in files0],
                 [os.path.join(tmp_path, f1) for f1 in files1],
                 repeat(os.path.join(data_dir, 'adaptive_gls.csv'), len(indices)),  # path to lookup table
                 repeat(tmp_path, len(indices))))  # directory for temporary output

with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(poolvcf.pysam_pooler_gt, args1))

print('\r\nTime elapsed --> ', timeit.default_timer() - start)

### CONCATENATE RESULTS
files2 = ['{}.gz'.format(f1) for f1 in files1]
pybcf.concat(files2, basin, tmp_path)
pybcf.sort(basin, tmp_path)
pybcf.bgzip(basin, basout, tmp_path)
pybcf.index(basout, tmp_path)
delete_file(basin)
shutil.copy(basout, os.path.expanduser(argsin.pathout))
print('\r\nTime elapsed --> ', timeit.default_timer() - start)
