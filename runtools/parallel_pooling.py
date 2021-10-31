import sys, os
import numpy as np
from itertools import starmap, repeat
import shutil
import multiprocessing as mp
import argparse
import timeit

home_dir = os.path.expanduser("~")
proj_dir = os.path.join(home_dir, '../')

# force PYTHONPATH to look into the project directory for modules
rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import poolvcf
from genotypooler.poolSNPs import pybcf

from genotypooler.persotools.files import delete_file, mkdir

'''
Parallelized file processing simulating pooling with adaptive genotype likelihoods.
For large VCF-files with some dozen of thousands or a million of variants.
If million of variants, it is recommened to run this script
on a cluster of servers e.g. UPPMAX resources (https://www.uppmax.uu.se/) .

Steps:
* Read main VCF file and write chunks: bash script based on bcftools
* Simulate pooling on chunks and decode into adaptive GP
* Merge pooled chunks back to read-for-use compressed VCF file

Usage (from genotypooler/runtools folder): 
$ python3 parallel_pooling.py ../data/IMP.chr20.snps.gt.vcf.gz ../data/IMP.chr20.pooled.snps.gl.vcf.gz 4
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
    tmp_path = os.path.join(os.path.dirname(argsin.pathin), 'tmp')  # os.path.join(proj_dir, '../data/tmp')
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

try:
    assert os.path.exists('snic')
    snic_proj = '/crex/proj/snic2019-8-216/private'
except AssertionError:
    snic_proj = os.path.dirname(argsin.pathin)  # os.path.join(proj_dir, '../data')
print('SNIC PROJ: {}'.format(snic_proj))

print('\n'.ljust(80, '*'))
print('Number of cpu to be used = {}'.format(nb_cores))
print('Input file = {}'.format(os.path.expanduser(argsin.pathin)))
print('Output file = {}'.format(os.path.expanduser(argsin.pathout)))
print('\n'.rjust(80, '*'))

### SIMULATE POOLING ON PACKED DATA CHUNKS
data_dir = snic_proj
os.chdir(data_dir)
fingz = argsin.pathin  # os.path.expanduser(argsin.pathin)  # /home/usr
basin = os.path.basename(fingz).rstrip('.gz').replace('.gt', '.gl')
basout = os.path.basename(argsin.pathout)


start = timeit.default_timer()
prename = 'pack'

files0 = [f for f in os.listdir(tmp_path)]
# keep only chunks to process and index the files
for f0 in files0:
    if not f0.startswith(prename) or f0.endswith('.csi'):
        files0.remove(f0)
    else:
        pybcf.index(f0, tmp_path)

indices = np.arange(len(files0))

print('\r\n{} files found will be pooled'.format(len(files0)).ljust(80, '.'))
print([os.path.join(tmp_path, f0) for f0 in files0])

files1 = ['pooled.{}'.format(f0).rstrip('.gz').replace('.gt', '.gl') for f0 in files0]

args1 = list(zip([f0 for f0 in files0],
                 [f1 for f1 in files1],
                 repeat(os.path.join(data_dir, 'adaptive_gls.csv'), len(indices)),  # path to lookup table
                 repeat(tmp_path, len(indices))))  # directory for temporary output

with mp.Pool(processes=nb_cores) as mpool:
    _ = all(mpool.starmap(poolvcf.pysam_pooler_gp, args1))

print('\r\nTime elapsed --> ', timeit.default_timer() - start)

### CONCATENATE RESULTS
files2 = ['{}.gz'.format(f1) for f1 in files1]
pybcf.concat(files2, basin, tmp_path)
pybcf.sort(basin, tmp_path)
pybcf.bgzip(basin, basout, tmp_path)
pybcf.index(basout, tmp_path)
os.chdir(data_dir)
delete_file(os.path.join(tmp_path, basin))
# shutil.copy(basout, os.path.join(tmp_path, basout))
print('\r\nTime elapsed --> ', timeit.default_timer() - start)
