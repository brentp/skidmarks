"""

Skid Marks: Check for runs in sequences
----------------------------------------

Q: how do you check for runs?
A: look for skidmarks.

This module implements some functions to check a sequence for randomness.
in some cases, it is assumed to be a binary sequence (not only 1's and 0's
but containing only 2 distinct values.
Any feedback or improvements are welcomed

    >>> from skidmarks import gap_test, wald_wolfowitz, auto_correlation, serial_test

"""

import math
from scipy.stats import zprob, linregress, chisquare
from itertools import groupby
# yoavram: Bug fix for Python 3 as suggested in https://github.com/nschloe/matplotlib2tikz/issues/20
try:
    from itertools import izip
except ImportError:
    izip = zip
import numpy as np
import collections

# yoavram: fix string import stuff in python 3, see https://github.com/oxplot/fysom/issues/1
try:
    unicode = unicode
except NameError:
    # 'unicode' is undefined, must be Python 3
    str = str
    unicode = str
    bytes = bytes
    basestring = (str,bytes)
else:
    # 'unicode' exists, must be Python 2
    str = str
    unicode = unicode
    bytes = str
    basestring = basestring
    

def wald_wolfowitz(sequence):
    """
    implements the wald-wolfowitz runs test:
    http://en.wikipedia.org/wiki/Wald-Wolfowitz_runs_test
    http://support.sas.com/kb/33/092.html

    :param sequence: any iterable with at most 2 values. e.g.
                     '1001001'
                     [1, 0, 1, 0, 1]
                     'abaaabbba'

    :rtype: a dict with keys of 
        `n_runs`: the number of runs in the sequence 
        `p`: the support to reject the null-hypothesis that the number of runs 
             supports a random sequence
        `z`: the z-score, used to calculate the p-value 
        `sd`, `mean`: the expected standard deviation, mean the number of runs, 
                      given the ratio of numbers of 1's/0's in the sequence

    >>> r = wald_wolfowitz('1000001')
    >>> r['n_runs'] # should be 3, because 1, 0, 1
    3

    >>> r['p'] < 0.05 # not < 0.05 evidence to reject Ho of random sequence
    False

    # this should show significance for non-randomness
    >>> li = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    >>> wald_wolfowitz(li)['p'] < 0.05
    True

    """
    R = n_runs = sum(1 for s in groupby(sequence, lambda a: a))

    n = float(sum(1 for s in sequence if s == sequence[0]))
    m = float(sum(1 for s in sequence if s != sequence[0]))

    # expected mean runs
    ER = ((2 * n * m ) / (n + m)) + 1
    # expected variance runs
    VR = (2 * n * m * (2 * n * m - n - m )) / ((n + m)**2 * (n + m - 1)) 
    O = (ER - 1) * (ER - 2) / (n + m - 1.)
    assert VR - O < 0.001, (VR, O)

    SD = math.sqrt(VR)
    # Z-score
    Z = (R - ER) / SD

    return {'z': Z, 'mean': ER, 'sd': SD, 'p': zprob(Z), 'n_runs': R}


def auto_correlation(sequence):
    """
    test for the autocorrelation of a sequence between t and t - 1
    as the 'auto_correlation' it is less likely that the sequence is
    generated randomly.
    :param sequence: any iterable with at most 2 values that can be turned
                     into an integer via int() . e.g.
                     '1001001'
                     [1, 0, 1, 0, 1]
    :rtype: returns a dict of the linear regression stats of sequence[1:] vs.
            sequence[:-1]

    >>> result = auto_correlation('00000001111111111100000000')
    >>> result['p'] < 0.05
    True
    >>> result['auto_correlation']
    0.83766233766233755

    """
    if isinstance(sequence, basestring):
        sequence = map(int, sequence)
    seq = np.array(sequence, dtype=np.int)
    dseq = np.column_stack((seq[1:], seq[:-1]))
    slope, intercept, r, ttp, see = linregress(seq[1:], seq[:-1])
    cc = np.corrcoef(dseq, rowvar=0)[0][1]
    return {'slope': slope, 'intercept': intercept, 'r-squared': r ** 2,
            'p': ttp, 'see': see, 'auto_correlation': cc}


def _avg_shannon_entropy(sequence):
    if isinstance(sequence, basestring):
        sequence = map(int, sequence)
    pm = [1 for s in sequence if s == sequence[0]]
    p = len(pm) / float(len(sequence))
    return sum([-p* math.log(p, 2) for _ in pm]) / (len(pm) or 1 / 2.0)



def serial_test(sequence):
    """
    serial test tests for randomness by looking at conversions from 1 digit to the next
    a low p-value in the returned dict indicates the strength that the null-hypothesis of a
    random sequence may be rejected.

    http://books.google.com/books?id=EIbxfCGfzgcC&lpg=PA141&ots=o-8ymmqbs9&pg=PA142#v=onepage&q=&f=false

    :param sequence: any iterable with at most 2 values that can be turned
                     into an integer via int() . e.g.
                     '1001001'
                     [1, 0, 1, 0, 1]
    :rtype: returns dict of {'chi': <chisquare value>, 'p': <p-value of said chisquare>}

    >>> serial_test('101010101111000')
    {'chi': 1.4285714285714286, 'p': 0.69885130769248427}

    >>> serial_test('110000000000000111111111111')
    {'chi': 18.615384615384617, 'p': 0.00032831021826061683}

    """
    #if isinstance(sequence, basestring): sequence = map(int, sequence)
    pairwise = izip(sequence[1:], sequence[:-1])
    d = collections.defaultdict(int)
    for k in pairwise: d[k] += 1
    # order doesnt matter because the expected are all the same.
    obs = np.array(list(d.values()))
    exp = np.ones_like(obs) * obs.mean()

    chi, pval =  chisquare(obs, exp)
    return {'chi': chi, 'p': pval}


def gap_test(sequence, item=None):
    """
    http://books.google.com/books?id=EIbxfCGfzgcC&lpg=PA141&ots=o-8ymmqbs9&pg=PA142#v=onepage&q=&f=false

    check for randomness in the sequence using distance (in the sequence) between `items`.

    the gap test for runs. takes a sequence and option `item`
    :param item: is used to test for gaps.
    :rtype: dict of pvalue, chi-square and the `item`
    >>> gap_test('100020001200000')
    {'chi': 756406.99909855379, 'item': '1', 'p': 0.0}

    >>> gap_test('101010111101000')
    {'chi': 11.684911193438811, 'item': '1', 'p': 0.23166089118674466}


    gap_test() will default to looking for gaps between the first value in
    the sequence (in this case '1') and each later occurrence. use the `item`
    kwarg to specify another value.

        >>> gap_test('101010111101000', item='0')
        {'chi': 11.028667632612191, 'item': '0', 'p': 0.27374903509732523}


    """

    if item is None: item=sequence[0] 

    seq = [1 if s == item else 0 for s in sequence]
    assert 0 < sum(seq) < len(sequence), \
            "must have some values of item: %s" % item

    observed_gaps = [len(list(li)) for g, li in groupby(seq, lambda a: a) if g != 1]
    ogaps = [observed_gaps.count(i) \
                     for i in range(max(10, max(observed_gaps)))]

    l = float(len(sequence))
    exp = float(sum(seq)) / l
    egaps = [l * (exp ** ii) for ii in range(1, len(ogaps) + 1)]
    chi, pval = chisquare(np.array(ogaps), np.array(egaps))
    return {'chi': chi, 'p': pval, 'item': item}

def test_suite():
    import doctest
    suite = doctest.DocFileSuite(__file__, module_relative=False)
    return suite
    
if __name__ == "__main__":

    import doctest
    doctest.testmod()
