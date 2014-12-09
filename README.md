Skid Marks: Check for runs in sequences
---------------------------------------

  Q: how do you check for runs?

  A: look for skidmarks.

This module implements some functions to check a sequence for randomness.
in some cases, it is assumed to be a binary sequence (not only 1's and 0's
but containing only 2 distinct values).
Any feedback, improvements, additions are welcomed.

    >>> from skidmarks import gap_test, wald_wolfowitz, auto_correlation, serial_test


Wald-Wolfowitz
--------------

http://en.wikipedia.org/wiki/Wald-Wolfowitz\_runs\_test

http://support.sas.com/kb/33/092.html

    >>> r = wald_wolfowitz('1000001')
    >>> r['n_runs'] # should be 3, because 1, 0, 1
    3

    >>> r['p'] < 0.05 # not < 0.05 evidence to reject Ho of random sequence
    False

# this should show significance for non-randomness
    >>> li = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    >>> wald_wolfowitz(li)['p'] < 0.05
    True



Autocorrelation
---------------

    >>> result = auto_correlation('00000001111111111100000000')
    >>> result['p'] < 0.05
    True

    >>> result['auto_correlation']
    0.83766233766233755


Serial Test
-----------

http://books.google.com/books?id=EIbxfCGfzgcC&lpg=PA141&ots=o-8ymmqbs9&pg=PA142#v=onepage&q=&f=false

    >>> serial_test('101010101111000')
    {'chi': 1.4285714285714286, 'p': 0.69885130769248427}

    >>> serial_test('110000000000000111111111111')
    {'chi': 18.615384615384617, 'p': 0.00032831021826061683}


Gap Test
--------

http://books.google.com/books?id=EIbxfCGfzgcC&lpg=PA141&ots=o-8ymmqbs9&pg=PA142#v=onepage&q=&f=false

    >>> gap_test('100020001200000')
    {'chi': 756406.99909855379, 'item': '1', 'p': 0.0}

    >>> gap_test('101010111101000')
    {'chi': 11.684911193438811, 'item': '1', 'p': 0.23166089118674466}

gap\_test() will default to looking for gaps between the first value in
the sequence (in this case '1') and each later occurrence. use the `item`
kwarg to specify another value.

    >>> gap_test('101010111101000', item='0')
    {'chi': 11.028667632612191, 'item': '0', 'p': 0.27374903509732523}
