import os
import numpy as np
import math
from scipy.linalg import block_diag
import itertools
import pandas as pd
from typing import *


"""
Classes and methods for simulating pooling on genotypes for chosen samples at one variant.
Can handle a single block or repeated blocks (Non Overlapping Repeated Block Design.
Samples are assumed to be sorted be consecutive pools by row-major order.
A variant is represented by a NumPy array with GT data for each sample.
"""


class Design(np.ndarray):
    """
    Design matrix and pooling design
    """

    def __new__(cls,
                shape: np.ndarray = np.asarray([4, 4]),
                id_len: int = 8,
                pools_nb: int = 8,
                pools_size: int = 4,
                blocks: int = 1) -> np.ndarray:
        """
        Define the basic structure for a pool i.e.
        a squared matrix to fill with the variables IDs/GT/GL.
        :param shape: tuple, shape of the pool
        :param id_len: max number of char of the variables IDs
        :param pools_nb: number of pools per block
        :param pools_size: pool's size within a block
        :param blocks: number of repeated blocks
        :return: a matrix with dims 'shape', filled with str types
        """
        cls.id_len = id_len
        #id = 'U' + str(cls.id_len)
        cls.pools_nb = pools_nb
        cls.pools_size = pools_size
        cls.blocks = blocks
        return np.empty_like(super(Design, cls).__new__(cls, shape), dtype=int)

    @property
    def matrix(self) -> np.ndarray:
        """
        That function is not intended to be called explicitly.
        :param random: bool for dispatching idv randomly in the matrix?
        :return: design matrix. Numpy array.
        """
        pools_size: int = self.pools_size
        m: np.ndarray = np.zeros((self.pools_nb, self.size), dtype=int)
        for i in range(int(self.pools_nb/self.ndim)):
            j = i * pools_size
            m[i, j:j+pools_size] = [1]*pools_size
        for i in range(int(self.pools_nb/self.ndim), self.pools_nb):
            j = i - pools_size
            m[i, [j+k*pools_size for k in range(pools_size)]] = 1
        b = itertools.repeat(m, self.blocks)
        M = block_diag(*b)
        return M

    def __array_finalize__(self, obj: object) -> None:
        """
        Constructor needed for subclassing NumPy arrays.
        See online documentation.
        :param obj:
        :return:
        """
        if obj is None: return
        self.info = getattr(obj, 'info', None)


class Encoder(object):
    """
    Simulate encoding step in pooling
    """
    def __init__(self, design: object, format: str = 'gt'):
        self.ds = design
        assert format == 'gt'
        self.fmt = format

    def encode(self, variant: np.ndarray):
        """
        :param variant: sum of alleles in GT format, binary array has shape (N,)
        where N is the number of samples
        """
        scores: np.ndarray = variant  # np.apply_along_axis(sum, axis=-1, arr=variant)
        pooled_gt = np.dot(self.ds,
                           np.transpose(scores)).reshape((1, self.ds.shape[0], 1))
        pooled_gt = np.broadcast_to(pooled_gt, (1, self.ds.shape[0], 2))
        p = np.apply_along_axis(self.gt_converter, axis=-1, arr=pooled_gt)

        return p  # list of gt for the 8 pools from design matrix

    def gt_converter(self, a: np.ndarray) -> np.ndarray:
        """
        Formats pooled scores into individual GT.
        :param a: score from matrix-vector pooling
        :return: pool's true genotype with phase
        """
        max_score = self.ds[0, :].sum() * 2  # diallelic markers assumed
        if np.all(a == 0):  # RR * RR * RR * RR
            gt = [0, 0]
        elif np.all(a == max_score):  # AA * AA * AA * AA
            gt = [1, 1]
        else:
            gt = [1, 0]
        return np.asarray(gt)


class DictBlockDecoder(object):
    """
        Simulate edeoding step in pooling.
        Proceed block-wise.
        This version is based on the use of a dictionary as lookup table with the adaptive GP values.
        """

    def __init__(self, design1block: np.ndarray, lookup: dict, format: str = 'gt'):
        self.ds = design1block  # block-wise
        if lookup is not None:
            self.dict_gl = lookup
        else:
            self.dict_gl = None
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
        Decode to adaptive log-GP if a lookup table with specific values is provided,
        else uniform random log-GP = (-0.47712, -0.47712, -0.47712) i.e. GP = (0.33333, 0.33333, 0.33333)
        """
        scores: np.ndarray = pooled.flatten()  # np.apply_along_axis(sum, axis=-1, arr=pooled).flatten()
        rowcolcounts = self.rowcolcounts(scores)
        colcross = np.apply_along_axis(np.multiply, 1, self.ds.T, scores)  # (16, 8)
        masks = np.ma.masked_where(self.ds.T, colcross)
        # Returns and sorts genotypes of pools cross section for an individual
        crosses = np.sort(colcross[masks.mask].reshape((self.ds.shape[1], 2)))  # (self.ds.shape[1], 2) = (16, 2) and 2 is the weight of the design
        if self.dict_gl is not None:
            unknown = [self.dict_gl[tuple([*rowcolcounts, *crs])] for crs in crosses]
        else:
            unknown = []
            for crs in crosses:
                if np.equal(crs, [0, 0]).all() or np.equal(crs, [0, 1]).all() or np.equal(crs, [1, 0]).all():
                    unknown.append([0., -12., -12.])
                elif np.equal(crs, [2, 2]).all() or np.equal(crs, [2, 1]).all() or np.equal(crs, [1, 2]).all():
                    unknown.append([-12., -12., 0.])
                else:
                    unknown.append([-0.47712, -0.47712, -0.47712])
        decoded_gp = np.asarray(unknown)

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
        #TODO: if GP and None lookup -> decode into 0.33, 0.33, 0.33


    def decode_genotypes_gt(self, pooled: np.ndarray) -> np.ndarray:
        """
        Recomputes true genotypes of samples with/without pooling/missing data
        :param pooled: sum of alleles of pooled true genotypes (unpohased)
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
        Decode to adaptive log-GP if a lookup table with specific values is provided,
        else uniform random log-GP = (-0.47712, -0.47712, -0.47712) i.e. GP = (0.33333, 0.33333, 0.33333)
        :return: individual samples genotypes (genotype likelihoods)
        """
        scores: np.ndarray = pooled.flatten()  # np.apply_along_axis(sum, axis=-1, arr=pooled).flatten()
        rowcolcounts = self.rowcolcounts(scores)
        colcross = np.apply_along_axis(np.multiply, 1, self.ds.T, scores)  # (16, 8)
        masks = np.ma.masked_where(self.ds.T, colcross)
        # Returns and sorts genotypes of pools cross section for an individual
        crosses = np.sort(colcross[masks.mask].reshape(
            (self.ds.shape[1], 2)))  # (self.ds.shape[1], 2) = (16, 2) and 2 is the weight of the design
        # this sorts colcross coordinates only (as sorted in the csv table too)
        keys = np.asarray([[*rowcolcounts, *crs] for crs in crosses])

        if self.V is not None and self.D is not None:
            unknown = np.apply_along_axis(self.multidecoder_gp, axis=-1, arr=keys)
        else:
            unknown = []
            for crs in crosses:
                if np.equal(crs, [0, 0]).all() or np.equal(crs, [0, 1]).all() or np.equal(crs, [1, 0]).all():
                    unknown.append([0., -12., -12.])
                elif np.equal(crs, [2, 2]).all() or np.equal(crs, [2, 1]).all() or np.equal(crs, [1, 2]).all():
                    unknown.append([-12., -12., 0.])
                else:
                    unknown.append([-0.47712, -0.47712, -0.47712])
        decoded_gl = np.asarray(unknown)

        return decoded_gl

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

    def multidecoder_gp(self, k: np.ndarray) -> np.ndarray:
        dkey = get_dummy_key(k)
        gidx = np.digitize(dkey.dot(self.D.transpose()), [len(k)])
        gp = gidx.dot(self.V)
        return gp


class Decoder(object):
    """
        Simulate deceoding step in pooling.
        Proceed block-wise.
        This version is based on the use of dot-products of NumPy arrays for
         representing the lookup table with the adaptive GP values.
        """

    def __init__(self, design_matrix: np.ndarray, lookup_keys: np.ndarray, lookup_vals: np.ndarray, format: str = 'gt'):
        self.dm = design_matrix  # matrix for all blocks
        self.ds1 = Design()  # single block
        self.dm1 = self.ds1.matrix
        self.fmt = format.upper()
        assert (self.fmt == 'GT' or self.fmt == 'GP'), 'Pooling to other formats than GT or GP not implemented'
        self.D, self.V = lookup_keys, lookup_vals
        #TODO: if GP and None lookup -> decode into 0.33, 0.33, 0.33


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
        Decode to adaptive log-GP if a lookup table with specific values is provided,
        else uniform random log-GP = (-0.47712, -0.47712, -0.47712) i.e. GP = (0.33333, 0.33333, 0.33333)
        """
        # Reshape for NumPy style np.apply_along_axis
        scores: np.ndarray = pooled.reshape((self.n_blocks, self.dm1.shape[0]))
        rowcolcounts = np.apply_along_axis(self.rowcolcounts, axis=-1, arr=scores)
        # Gets scores of the pair of pools intersecting at each sample
        colcross = np.apply_along_axis(np.multiply, 1, self.dm.T, pooled).squeeze()
        # Masks non-zeros i.e. pair of pools scores
        masks = np.ma.masked_where(self.dm.T, colcross)
        # Returns and sorts genotypes of pools cross section for an individual
        crosses = np.sort(colcross[masks.mask].reshape((self.dm.shape[1], 2)))  # 2 is the weight of the design
        # this sorts colcross coordinates only (as sorted in the csv table too)
        poolscounts = np.tile(rowcolcounts,
                              self.dm1.shape[1]).reshape(self.dm.shape[1], rowcolcounts.shape[1])
        keys = np.concatenate([poolscounts, crosses], axis=1)
        if self.V is not None and self.D is not None:
            decoded_gp = np.apply_along_axis(self.multidecoder_gp, axis=-1, arr=keys)
        else:
            unknown = []
            for crs in crosses:
                if np.equal(crs, [0, 0]).all() or np.equal(crs, [0, 1]).all() or np.equal(crs, [1, 0]).all():
                    unknown.append([0., -12., -12.])
                elif np.equal(crs, [2, 2]).all() or np.equal(crs, [2, 1]).all() or np.equal(crs, [1, 2]).all():
                    unknown.append([-12., -12., 0.])
                else:
                    unknown.append([-0.47712, -0.47712, -0.47712])
            decoded_gp = np.asarray(unknown)

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

    def multidecoder_gp(self, k: np.ndarray) -> np.ndarray:
        dkey = get_dummy_key(k)
        gidx = np.digitize(dkey.dot(self.D.transpose()), [len(k)])
        gp = gidx.dot(self.V)
        return gp


def load_lookup_table(path: str) -> pd.DataFrame:
    """
    Return table of adaptive GP values as a DataFrame.
    """
    if path is None:
        df = None
    else:
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

    D = D.values
    V = V.values
    return D, V


def get_dummy_key(k) -> np.ndarray:
    """
    Converts the key array to a categorical "dummy" array
    """
    ds = Design()
    rg = ds.pools_size + 1
    strides = np.asarray([rg, rg, rg, rg, rg, rg, 3, 3])  # 3 = genotypes RR, RA/AR, AA
    dumkey = np.zeros((strides.sum(),), dtype=int)
    idx = 0
    for i in range(len(k)):
        dumkey[idx + k[i]] = 1
        idx = idx + strides[i]
    return dumkey


def load_lookup_dict(path: str, log10: bool = True) -> dict:
    """
    Read adaptive GP from a csv file and return them as log-GP in a dictionary.
    Dictionary keys combine the pools' genotypes in a block and pairs of intersecting pools.
    """
    if path is None:
        df2dict = None
    else:
        df = pd.read_csv(os.path.abspath(path),
                         header=None,
                         names=['rowsrr', 'rowsra', 'rowsaa', 'colsrr', 'colsra', 'colsaa',
                                'n', 'm',
                                'rr', 'ra', 'aa']
                         )
        # log10func = lambda x: math.log10(x) if x > 1e-05 else -5.0
        log10func = lambda x: math.log10(x) if x > 1e-12 else -12.0
        if log10:
            df['rr'] = df['rr'].apply(log10func)
            df['ra'] = df['ra'].apply(log10func)
            df['aa'] = df['aa'].apply(log10func)
        df2dict = dict(((int(rwrr), int(rwra), int(rwaa), int(clrr), int(clra), int(claa),
                         int(n), int(m)),
                        [rr, ra, aa]) for rwrr, rwra, rwaa, clrr, clra, claa,
                                          n, m,
                                          rr, ra, aa in df.itertuples(index=False, name=None))
    return df2dict


def blocks_decoder(nB, v, step, lookup_keys, lookup_vals, dec_fmt: str):
    """
    Decodes nB blocks from a NORB pooling design
    :param step: number of pools per block
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
    Decodes nB blocks from a NORB pooling design
    :param step: number of samples per block
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

