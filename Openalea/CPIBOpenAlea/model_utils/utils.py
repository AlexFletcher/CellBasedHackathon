
"""This module contains a small number of general utility functions"""

def pairs(lst):
    """Generaterator which returns all the sequential pairs of lst,
    treating lst as a circular list
    
    >>> list(pairs([1,2,3]))
    [(1, 2), (2, 3), (3, 1)]
    """     

    i = iter(lst)
    first = prev = i.next()
    for item in i:
        yield prev, item
        prev = item
    yield item, first

def find_pred(s, pred):
    """ Returns index of first item in the iterable s for which pred(s) is true
    
    >>> find_pred([1,2,3], lambda x: x==2)
    1
    >>> find_pred([1,2,3], lambda x: x==4)
    """
    return next((i for i, v in enumerate(s) if pred(v)), None)


def sgn(x):
    """ Sign of a float
    
    >>> sgn(2.3)
    1
    >>> sgn(0.0)
    0
    >>> sgn(-3.7)
    -1
    """
    return cmp(x, 0.0)

def sort_pair(p):
    """ Sorts a pair of values (returning a tuple)
    
    >>> sort_pair((1,2))
    (1,2)
    >>> sort_pair((2,1))
    (1,2)
    """
    return p if (p[1]>p[0]) else (p[1], p[0])


