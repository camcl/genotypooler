import os
import numpy as np

from genotypooler.poolSNPs.pooler import Design, Encoder

import pandas as pd
from numba import jit, jitclass
from numba import int32, float32, float64, types
import timeit
from typing import *


"""
Classes and methods for simulating pooling on genotypes for chosen samples at one variant.
Can handle a single block or repeated blocks (Non Overlapping Repeated Block Design.
Samples are assumed to be sorted be consecutive pools by row-major order.
A variant is represented by a NumPy array with GT data for each sample.
"""

@jit(nopython=True, forceobj=True, locals={'countrr': int32, 'countra': int32, 'countaa': int32})
def rowcolcounter(a: np.ndarray) -> np.ndarray:
    """
    Count number of pooled RR|RA|AA genotypes over all rows and columns of a pooled matrix
    :param a: score i.e. trinary-encoded true genotypes
    :return: counts of genotypes for the rows and columns
    """
    # a should be scores
    count_rr = np.where(a == 0, 1, 0)
    count_aa = np.where(a == 2, 1, 0)
    count_ra = np.where(a == 1, 1, 0)

    rowcolcounts = np.zeros((3*2,), dtype=int)  # counts number of RR, RA, AA rows, same for columns
    rowcolcounts[0] = np.sum(count_rr[:4])
    rowcolcounts[3] = np.sum(count_rr[4:])
    rowcolcounts[1] = np.sum(count_ra[:4])
    rowcolcounts[4] = np.sum(count_ra[4:])
    rowcolcounts[2] = np.sum(count_aa[:4])
    rowcolcounts[5] = np.sum(count_aa[4:])

    return rowcolcounts


@jit(nopython=True, forceobj=True, locals={'gidx': float64, 'gp': float64})
def multidecoder_gp(dkey: np.ndarray, lookkeys: np.ndarray, lookvals: np.ndarray) -> np.ndarray:
    """gets GP from lookup table as arrays"""
    # gidx = np.floor_divide(dkey.dot(lookkeys.T), np.sum(dkey))  # np.sum(dkey) = 8
    # gidx = np.digitize(dkey.dot(lookkeys.T), [np.sum(dkey)])
    gidx = dkey.dot(lookkeys.T) // np.sum(dkey)
    gp = gidx.dot(lookvals)
    return gp


def get_dict_val(k, d: dict):
    return d[tuple(k)]


class DictBlockDecoder(object):
    """
        Simulate edeoding step in pooling.
        Proceed block-wise.
        This version is based on the use of a dictionary as lookup table with the adaptive GP values.
        """

    def __init__(self, design1block: np.ndarray, lookup: dict, format: str = 'gt'):
        self.ds = design1block  # block-wise
        self.dict_gl = lookup
        self.fmt = format.upper()
        assert (self.fmt == 'GT' or self.fmt == 'GP'), 'Pooling to other formats than GT or GP not implemented'

    def decode_genotypes_gt(self, pooled: np.ndarray) -> np.ndarray:
        """
        Recomputes true genotypes of samples with/without pooling/missing data
        :param pooled: sum of alleles of pooled true genotypes (unpohased)
        :return: individual samples genotypes (true genotype unphased)
        """
        count_alt: np.ndarray = np.where(pooled >= 1, 1, 0)
        count_ref: np.ndarray = np.where(pooled <= 1, 1, 0)
        # 2 = value max of allele dosage with diploid diallelic markers

        alt_row: int = np.sum(count_alt[4:])
        alt_col: int = np.sum(count_alt[4:])
        ref_row: int = np.sum(count_ref[:4])
        ref_col: int = np.sum(count_ref[4:])

        nb_alt: int = alt_row + alt_col
        nb_ref: int = ref_row + ref_col

        encoded = np.dot(pooled,
                         self.ds).reshape(1, self.ds.shape[1], 1)
        b = np.broadcast_to(encoded, (1, self.ds.shape[1], 2))
        if nb_alt == 0:
            decoded_gt = np.zeros_like(b)
        elif nb_ref == 0:
            aa = np.array([1, 1])
            decoded_gt = np.tile(aa, self.ds.shape[1]).reshape((1, self.ds.shape[1], 2))
        elif nb_alt == 2:
            decoder: Callable = lambda x: [1, -1] if np.all(x == 2) else [0, 0]
            # np.all() because of b.shape
            decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
        elif nb_ref == 2:  # symmetric case with ref allele in only 2 pools: individual is RR or RA
            decoder: Callable = lambda x: [0, -1] if np.all(x == 2) else [1, 1]
            decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
        else:  # nb_alt > 2 and nb_alt < 8: # nb_alt = 2*n with n > 1
            decoded_gt = np.apply_along_axis(self.multidecoder_gt, axis=-1, arr=b)

        return decoded_gt

    def decode_genotypes_gp(self, pooled: np.ndarray) -> np.ndarray:
        """

        """
        scores: np.ndarray = pooled.flatten()  # np.apply_along_axis(sum, axis=-1, arr=pooled).flatten()
        rowcolcounts = rowcolcounter(scores)
        colcross = np.apply_along_axis(np.multiply, 1, self.ds.T, scores)  # (16, 8)
        masks = np.ma.masked_where(self.ds.T, colcross)
        # Returns and sorts genotypes of pools cross section for an individual
        crosses = np.sort(colcross[masks.mask].reshape((self.ds.shape[1], 2)))  # (self.ds.shape[1], 2) = (16, 2) and 2 is the weight of the design
        poolscounts = np.tile(rowcolcounts,
                              self.ds.shape[1]).reshape((self.ds.shape[1], rowcolcounts.shape[0]))
        kys = np.concatenate([poolscounts, crosses], axis=1)
        unknown = [self.dict_gl[tuple([*rowcolcounts, *crs])] for crs in crosses]
        decoded_gp = np.asarray(unknown)
        # decoded_gp = np.apply_along_axis(get_dict_val, -1, kys, self.dict_gl)

        return decoded_gp

    @staticmethod
    def rowcolcounts(a: np.ndarray) -> np.ndarray:
        """
        Count number of pooled RR|RA|AA genotypes over all rows and columns of a pooled matrix
        :param a: score i.e. trinary-encoded true genotypes
        :return: counts of genotypes for the rows and columns
        """
        # a should be scores
        count_rr: np.ndarray = np.where(a == 0, 1, 0)
        count_aa: np.ndarray = np.where(a == 2, 1, 0)
        count_ra: np.ndarray = np.where(a == 1, 1, 0)

        rowcolcounts: np.ndarray = np.zeros((3 * 2,), dtype=int)  # counts number of RR, RA, AA rows, same for columns
        rowcolcounts[0] = np.sum(count_rr[:4])
        rowcolcounts[3] = np.sum(count_rr[4:])
        rowcolcounts[1] = np.sum(count_ra[:4])
        rowcolcounts[4] = np.sum(count_ra[4:])
        rowcolcounts[2] = np.sum(count_aa[:4])
        rowcolcounts[5] = np.sum(count_aa[4:])

        return rowcolcounts

    @staticmethod
    def multidecoder_gt(a: np.ndarray) -> np.ndarray:
        """
        Decodes pooled scores into individual GT.
        :param a: score
        :return: true genotype with phase
        """
        if np.all(a == 2):  # RA * RA
            gt = [-1, -1]
        elif np.all(a == 1) or np.all(a == 0):  # RA * RR or RR * RR
            gt = [0, 0]
        else:
            gt = [1, 1]
        return np.asarray(gt)


class SingleBlockDecoder(object):
    """
        Simulate edeoding step in pooling.
        Proceed block-wise.
        This version is based on the use of dot-products of NumPy arrays for
         representing the lookup table with the adaptive GP values.
        """

    def __init__(self, design1block: np.ndarray, lookup_keys: np.ndarray, lookup_vals: np.ndarray, format: str = 'gt'):
        self.ds = design1block  # block-wise
        self.fmt = format.upper()
        assert (self.fmt == 'GT' or self.fmt == 'GP'), 'Pooling to other formats than GT or GP not implemented'
        self.D, self.V = lookup_keys, lookup_vals

    def decode_genotypes_gt(self, pooled: np.ndarray) -> np.ndarray:
        """
        Recomputes true genotypes of samples with/without pooling/missing data
        :param pooled: sum of alleles of pooled true genotypes (unphased)
        :return: individual samples genotypes (true genotype unphased)
        """
        count_alt: np.ndarray = np.where(pooled >= 1, 1, 0)
        count_ref: np.ndarray = np.where(pooled <= 1, 1, 0)

        alt_row: int = np.sum(count_alt[4:])
        alt_col: int = np.sum(count_alt[4:])
        ref_row: int = np.sum(count_ref[:4])
        ref_col: int = np.sum(count_ref[4:])

        nb_alt: int = alt_row + alt_col
        nb_ref: int = ref_row + ref_col

        encoded = np.dot(pooled,
                         self.ds).reshape(1, self.ds.shape[1], 1)
        b = np.broadcast_to(encoded, (1, self.ds.shape[1], 2))
        if nb_alt == 0:
            decoded_gt = np.zeros_like(b)
        elif nb_ref == 0:
            aa = np.array([1, 1])
            decoded_gt = np.tile(aa, self.ds.shape[1]).reshape((1, self.ds.shape[1], 2))
        elif nb_alt == 2:
            decoder: Callable = lambda x: [1, -1] if np.all(x == 2) else [0, 0]
            # np.all() because of b.shape
            decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
        elif nb_ref == 2:  # symmetric case with ref allele in only 2 pools: individual is RR or RA
            decoder: Callable = lambda x: [0, -1] if np.all(x == 2) else [1, 1]
            decoded_gt = np.apply_along_axis(decoder, axis=-1, arr=b)
        else:  # nb_alt > 2 and nb_alt < 8: # nb_alt = 2*n with n > 1
            decoded_gt = np.apply_along_axis(self.multidecoder_gt, axis=-1, arr=b)

        return decoded_gt

    def decode_genotypes_gp(self, pooled: np.ndarray) -> np.ndarray:
        """
        Recomputes genotypes likelihoods of samples with/without pooling/missing data
        :param samples_gl: samples' true genotypes with phase
        :param dict_gl: likelihoods values to set when encountering missing genotypes
        :return: individual samples genotypes (genotype likelihoods)
        """
        scores: np.ndarray = pooled.flatten()  # np.apply_along_axis(sum, axis=-1, arr=pooled).flatten()
        rowcolcounts = rowcolcounter(scores)
        colcross = np.apply_along_axis(np.multiply, 1, self.ds.T, scores)  # (16, 8)
        masks = np.ma.masked_where(self.ds.T, colcross)
        # Returns and sorts genotypes of pools cross section for an individual
        crosses = np.sort(colcross[masks.mask].reshape((self.ds.shape[1], 2)))  # (self.ds.shape[1], 2) = (16, 2) and 2 is the weight of the design
        # this sorts colcross coordinates only (as sorted in the csv table too)
        poolscounts = np.tile(rowcolcounts,
                              self.ds.shape[1]).reshape((self.ds.shape[1], rowcolcounts.shape[0]))
        kys = np.concatenate([poolscounts, crosses], axis=1)
        dkeys = np.apply_along_axis(get_dummy_key, -1, kys)
        unknown = np.apply_along_axis(multidecoder_gp, -1, dkeys, self.D, self.V)
        decoded_gl = np.asarray(unknown)

        return decoded_gl

    @staticmethod
    def multidecoder_gt(a: np.ndarray) -> np.ndarray:
        """
        Decodes pooled scores into individual GT.
        :param a: score
        :return: true genotype with phase
        """
        if np.all(a == 2):  # RA * RA
            gt = [-1, -1]
        elif np.all(a == 1) or np.all(a == 0):  # RA * RR or RR * RR
            gt = [0, 0]
        else:
            gt = [1, 1]
        return np.asarray(gt)


class Decoder(object):
    """
        Simulate edeoding step in pooling.
        Proceed block-wise.
        This version is based on the use of dot-products of NumPy arrays for
         representing the lookup table with the adaptive GP values.
        """

    # def __init__(self, design_matrix: np.ndarray, lookup_keys: np.ndarray, lookup_vals: np.ndarray, format: str = 'gt'):
    #     self.dm = design_matrix  # matrix for all blocks
    #     self.ds1 = Design()  # single block
    #     self.dm1 = self.ds1.matrix
    #     self.fmt = format.upper()
    #     assert (self.fmt == 'GT' or self.fmt == 'GP'), 'Pooling to other formats than GT or GP not implemented'
    #     self.D, self.V = lookup_keys, lookup_vals
    #     #TODO: if GP and None lookup -> decode into 0.33, 0.33, 0.33

    def __init__(self, design_matrix: np.ndarray, lookup: dict, format: str = 'gt'):
        self.dm = design_matrix  # matrix for all blocks
        self.ds1 = Design()  # single block
        self.dm1 = self.ds1.matrix
        self.dict_gl = lookup
        self.fmt = format.upper()
        assert (self.fmt == 'GT' or self.fmt == 'GP'), 'Pooling to other formats than GT or GP not implemented'


    @property
    def n_blocks(self):
        return self.dm.shape[1] // self.dm1.shape[1]

    def decode_genotypes_gt(self, pooled: np.ndarray) -> np.ndarray:
        """
        Recomputes true genotypes of samples with/without pooling/missing data
        :param pooled: sum of alleles of pooled true genotypes (unpohased)
        :return: individual samples genotypes (true genotype unphased)
        """
        where_alt:  np.ndarray = np.where(pooled >= 1, 1, 0)
        where_ref: np.ndarray = np.where(pooled <= 1, 1, 0)

        count_alt: np.ndarray = where_alt.reshape((1, self.n_blocks, self.dm1.shape[0])).sum(axis=-1)  # 1 count per block
        count_ref: np.ndarray = where_ref.reshape((1, self.n_blocks, self.dm1.shape[0])).sum(axis=-1)

        nb_alt: np.ndarray = np.repeat(count_alt, self.dm1.shape[1])  # 1 count per sample
        nb_ref: np.ndarray = np.repeat(count_ref, self.dm1.shape[1])

        encoded = np.dot(pooled, self.dm).reshape(1, self.dm.shape[1], 1)
        scores = np.full((5, self.dm.shape[1]), -1, dtype=int)
        # 5 rows = 1 row for nb_alt, 1 row for nb_ref, 1 row for encoded, 2 rows for the decoded genotype
        scores[0] = nb_alt
        scores[1] = nb_ref
        scores[2] = encoded.squeeze()
        decoded_gt = np.apply_along_axis(self.multidecoder_gt, axis=0, arr=scores)[3:]

        return decoded_gt.T

    @staticmethod
    def multidecoder_gt(a: np.ndarray) -> np.ndarray:
        """
        Decodes pooled scores into individual GT.
        :param a: score as quintet
        :return:
        """
        if a[0] == 0:  # all RR
            a[3] = 0
            a[4] = 0
        elif a[1] == 0:  # all AA
            a[3] = 1
            a[4] = 1
        elif a[0] == 2:   # all RR but 1 AA or RA
            if a[2] == 2:
                a[3] = 1
            else:
                a[3] = 0
                a[4] = 0
        elif a[1] == 2:   # symmetric case: all AA but 1 RR or RA
            if a[2] == 2:
                a[3] = 0
            else:
                a[3] = 1
                a[4] = 1
        else:  # mix of RR, AA, RA/AR
            if a[2] > 2:
                a[3] = 1
                a[4] = 1
            elif a[2] < 2:
                a[3] = 0
                a[4] = 0
        return a

    def decode_genotypes_gp(self, pooled: np.ndarray) -> np.ndarray:
        """

        """
        # Reshape for NumPy style np.apply_along_axis
        scores: np.ndarray = pooled.reshape((self.n_blocks, self.dm1.shape[0]))
        rowcolcounts = np.apply_along_axis(rowcolcounter, axis=-1, arr=scores)
        # Gets scores of the pair of pools intersecting at each sample
        colcross = np.apply_along_axis(np.multiply, 1, self.dm.T, pooled).squeeze()
        # Masks non-zeros i.e. pair of pools scores
        masks = np.ma.masked_where(self.dm.T, colcross)
        # Returns and sorts genotypes of pools cross section for an individual
        crosses = np.sort(colcross[masks.mask].reshape((self.dm.shape[1], 2)))  # 2 is the weight of the design
        # this sorts colcross coordinates only (as sorted in the csv table too)
        poolscounts = np.tile(rowcolcounts,
                              self.dm1.shape[1]).reshape(self.dm.shape[1], rowcolcounts.shape[1])
        kys = np.concatenate([poolscounts, crosses], axis=1)
        # dkeys = np.apply_along_axis(get_dummy_key, -1, kys)
        # decoded_gp = np.apply_along_axis(multidecoder_gp, -1, dkeys, self.D, self.V)
        # dict style
        decoded_gp = np.apply_along_axis(get_dict_val, -1, kys, self.dict_gl)

        return decoded_gp

    @staticmethod
    def rowcolcounts(a: np.ndarray) -> np.ndarray:
        """
        Count number of pooled RR|RA|AA genotypes over all rows and columns of a pooled matrix
        :param a: score i.e. trinary-encoded true genotypes
        :return: counts of genotypes for the rows and columns
        """
        # a should be scores
        count_rr: np.ndarray = np.where(a == 0, 1, 0)
        count_aa: np.ndarray = np.where(a == 2, 1, 0)
        count_ra: np.ndarray = np.where(a == 1, 1, 0)

        rowcolcounts: np.ndarray = np.zeros((3*2,), dtype=int)  # counts number of RR, RA, AA rows, same for columns
        rowcolcounts[0] = np.sum(count_rr[:4])
        rowcolcounts[3] = np.sum(count_rr[4:])
        rowcolcounts[1] = np.sum(count_ra[:4])
        rowcolcounts[4] = np.sum(count_ra[4:])
        rowcolcounts[2] = np.sum(count_aa[:4])
        rowcolcounts[5] = np.sum(count_aa[4:])

        return rowcolcounts

    # def multidecoder_gp(self, k: np.ndarray) -> np.ndarray:
    #     dkey = get_dummy_key(k)
    #     gidx = np.digitize(dkey.dot(self.D.transpose()), [len(k)])
    #     gp = gidx.dot(self.V)
    #     return gp


def load_lookup_table(path: str) -> pd.DataFrame:
    """
    Provides adaptive GL values as a DataFrame
    """
    df = pd.read_csv(os.path.join(path, 'adaptive_gls.csv'),
                     header=None,
                     names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                            'n', 'm',
                            'rr', 'ra', 'aa']
                     )
    return df


def get_lookup_arrays(path: str) -> Tuple[Type[np.ndarray], Any]:
    """
    Converts a lookup table to:
    * D: a categorical array encoding key as categorical ("dummy" array)
    * V: a value array where the position of any value matches its key position in D
    """
    T = load_lookup_table(path)
    T.drop_duplicates(inplace=True)
    T.reset_index(drop=True, inplace=True)

    dumlist = []
    V = T[T.columns[-3:]]
    for col in T.columns[:-3]:
        dumlist.append(pd.get_dummies(T[col], prefix=col))
    D = dumlist.pop(0)
    for df in dumlist:
        D = D.join(df)

    D = D.values.astype(float)
    V = V.values.astype(float)
    return D, V

@jit(nopython=True)
def get_dummy_key(k) -> np.ndarray:
    """
    Converts the key array to a categorical "dummy" array
    """
    rg = int(np.sum(k[:3]) + 1)  # 4 + 1 ex. sum(2, 2, 0) + 1
    strides = np.asarray([rg, rg, rg, rg, rg, rg, 3, 3])  # 3 = genotypes RR, RA/AR, AA
    dumkey = np.zeros((int(np.sum(strides)),))  #, dtype=int)
    idx = 0
    for i in range(len(k)):
        dumkey[idx + k[i]] = 1
        idx = idx + strides[i]
    return dumkey


def load_lookup_dict(path: str) -> dict:
    df = pd.read_csv(os.path.join(path, 'adaptive_gls.csv'),
                     header=None,
                     names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                            'n', 'm',
                            'rr', 'ra', 'aa']
                     )
    df2dict = dict(((int(rwrr), int(rwra), int(rwaa), int(clrr), int(clra), int(claa),
                     int(n), int(m)),
                    [rr, ra, aa]) for rwrr, rwra, rwaa, clrr, clra, claa,
                                      n, m,
                                      rr, ra, aa in df.itertuples(index=False, name=None))
    return df2dict


def blocks_decoder(nB, v, step, lookup_keys, lookup_vals, dec_fmt: str):
    """
    Decodes a single block from a NORB pooling design
    """
    ds = Design()
    dm = ds.matrix
    decoder = SingleBlockDecoder(dm, lookup_keys, lookup_vals, format=dec_fmt)
    res = []
    for b in range(nB):
        p = v.squeeze()[step*b:step*b + step]
        if decoder.fmt == 'GP':
            q = decoder.decode_genotypes_gp(p)
        elif decoder.fmt == 'GT':
            q = decoder.decode_genotypes_gt(p)
        res.append(q)
    return np.asarray(res).squeeze()


def dict_blocks_decoder(nB, v, step, lookup: dict, dec_fmt: str):
    """
    Decodes a single block from a NORB pooling design
    """
    ds = Design()
    dm = ds.matrix
    decoder = DictBlockDecoder(dm, lookup, format=dec_fmt)
    res = []
    for b in range(nB):
        p = v.squeeze()[step*b:step*b + step]
        if decoder.fmt == 'GP':
            q = decoder.decode_genotypes_gp(p)
        elif decoder.fmt == 'GT':
            q = decoder.decode_genotypes_gt(p)
        res.append(q)
    return np.asarray(res).squeeze()


def single_block_decoder(p, lookup_keys, lookup_vals, dec_fmt: str):
    ds = Design()
    dm = ds.matrix
    decoder = SingleBlockDecoder(dm, lookup_keys, lookup_vals, format=dec_fmt)
    if decoder.fmt == 'GP':
        q = decoder.decode_genotypes_gp(p)
    elif decoder.fmt == 'GT':
        q = decoder.decode_genotypes_gt(p)
    return q
