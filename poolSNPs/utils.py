from scipy.stats import *
import numpy as np
import pandas as pd

from genotypooler.persotools.debugging import *

dbg = MyPrintClass(True)

"""
Utils for data sets processing
NumPy or pandas DataFrame data sets
"""


def pgcd(a, b):
    """Computes the biggest common divider"""
    while b != 0:
        r = a % b
        a, b = b, r
    return a


def ppcm(a, b):
    """Computes the smallest common multiple"""
    if (a == 0) or (b == 0):
        return 0
    else:
        return (a*b)//pgcd(a, b)


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    else:
        return v / norm


#TODO: PEP8 refactoring
def sort_datasets(dflist: list, groups: list, df_aaf: pd.DataFrame) -> list:
    out = []
    df_aaf.reset_index(drop=True, inplace=True)
    for dfset in dflist:
        # Sort samples by population
        dfset.sort_index(axis=1, inplace=True)
        dfset.columns = groups
        dfset.sort_index(level='Population', axis=1, inplace=True)
        # Sort variants by AAF
        dfset.sort_index(axis=0, inplace=True) # sort by id
        dfset.reset_index(drop=True, inplace=True)
        dfset['af_info'] = df_aaf['af_info']
        # dfset['aaf_bin'] = df_aaf['aaf_bin']
        dfset.set_index([df_aaf['id'], 'af_info', 'aaf_bin'], drop=True, append=False, inplace=True)  # replace idx with multiidx (id sorted)
        # dfset.sort_index(level=['aaf_bin', 'af_info'], axis=0, inplace=True)
        dfset.sort_index(level=['af_info'], axis=0, inplace=True)
        dfset.reset_index(drop=True, inplace=True)
        out.append(dfset)
    return out


"""
Iterators vs Iterables
an iterable is an object that can return an iterator

    Examples: lists, strings, dictionaries, file connections

    An object with an associated iter() method

    Applying iter() to an iterable creates an iterator

an iterator is an object that keeps state and produces the next value when you call next() on it.

    Produces next value with next()
"""


def is_iterator(obj):
    if (hasattr(obj, '__iter__') and
            hasattr(obj, '__next__') and
            callable(obj.__iter__) and
            obj.__iter__() is obj):
        return True
    else:
        return False


def is_iterable(obj):
    if hasattr(obj, '__iter__') and callable(obj.__iter__):
        return True
    else:
        return False


def chunk(itr, chunk_size):
    """Generate sequences of `chunk_size` elements from `iterable`."""
    iterator = iter(itr)
    while True:
        chunk = []
        try:
            for _ in range(chunk_size):
                chunk.append(iterator.next())
            yield chunk
        except StopIteration:
            if chunk:
                yield chunk
            break

