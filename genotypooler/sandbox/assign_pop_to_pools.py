import os, sys
import numpy as np
from scipy.linalg import block_diag
import pandas as pd

rootdir = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, rootdir)

from genotypooler.poolSNPs import pooler


pop_size = 384
nb_blocks = 6
nb_pools_per_block = 8 + 8
nb_idv_per_pool = 8

# Seed pseudo-random generator
rng = np.random.default_rng(2023)

# Read list of individual IDs and shuffle it (possibly seed)
inlist = np.arange(1, pop_size + 1)
rng.shuffle(inlist)

# Define the design matrix for 8x8 NORB pooling design and 6 blocks (rows=pools, columns=samples)
design8x8 = pooler.Design(shape=np.asarray([nb_idv_per_pool, nb_idv_per_pool]),
                          id_len=3,
                          pools_nb=nb_pools_per_block,
                          pools_size=nb_idv_per_pool,
                          blocks=nb_blocks
                          )

# Construct data frame
df8x8 = pd.DataFrame(data=design8x8.matrix,
                     index=np.arange(1, nb_pools_per_block * nb_blocks + 1),
                     columns=inlist.astype(str))

outlist = df8x8.mul(inlist, axis=1).values
outdata = outlist[outlist != 0].reshape(nb_pools_per_block*nb_blocks,
                                        nb_idv_per_pool)

for k in range(nb_blocks):
    idv_in_block_k = np.unique(outdata[k*nb_pools_per_block:(k+1)*nb_pools_per_block,:nb_idv_per_pool])
    print(idv_in_block_k)


a = np.arange(0, nb_pools_per_block * nb_blocks)
# print(np.repeat(np.arange(1, nb_blocks + 1), nb_idv_per_pool ** 2))
block_pool_index = pd.MultiIndex.from_arrays([a // 16 + 1, a + 1], names=('block', 'pool'))
dfout = pd.DataFrame(data=outdata,
                     index=block_pool_index,
                     columns=list([f'sample{i}' for i in range(1, 8+1)])
                     ).reset_index()

dfout.to_csv('/home/camille/Desktop/pools_broccoli_avalo.csv', index=False)

dfsingles = pd.DataFrame(data=zip(np.repeat(np.arange(1, nb_blocks + 1), nb_idv_per_pool ** 2),
                                  np.arange(nb_pools_per_block * nb_blocks + 1,
                                            nb_pools_per_block * nb_blocks + pop_size + 1),
                                  inlist),
                         columns=['block', 'pool', 'sample1']
                         )

dfall = pd.concat([dfout, dfsingles], axis=0)\
    .fillna(0)\
    .astype(int)\
    .sort_values(by=['block', 'pool'])\
    .astype(str)\
    .replace('0', '')
print(dfall)

dfall.to_csv('/home/camille/Desktop/singles_and_pools_broccoli_avalo.csv',
             index=False,
             float_format='.0f'
             )
