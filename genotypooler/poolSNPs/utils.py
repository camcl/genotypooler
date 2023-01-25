import numpy as np

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

