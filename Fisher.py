import scipy.stats as ss
from scipy.special import gammaln as lgamma
from statlib.fexact import fisher_exact as f
from statlib.asa159 import rcont2
from statlib.asa205 import enum as rcont
import numpy as np
import logging
import os
import random


def fisher_exact(table, alternative="two-sided", hybrid=False, midP=False,
                 simulate_pval=False, replicate=2000, workspace=300,
                 attempt=2, seed=None):
    """Performs a Fisher exact test on a 2x2 contingency table.
    Parameters
    ----------
    table : array_like of ints
        A 2x2 contingency table.  Elements should be non-negative integers.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Which alternative hypothesis to the null hypothesis the test uses.
        Default is 'two-sided'.  Only used in the 2 x 2 case (with the scipy
        function). In every other case, the two-sided pval is returned.
    mult : int
        Specify the size of the workspace used in the network algorithm.
        Only used for non-simulated p-values larger than 2 x 2 table.
        You might want to increase this if the p-value failed!
    hybrid : bool
        Only used for larger than 2 x 2 tables, in which cases it indicates 
        whether the exact probabilities (default) or a hybrid approximation 
        thereof should be computed.
    midP : bool
        Use this to enable mid-P correction. Could lead to slow computation.
        This is not applicable for simulation p-values. `alternative` cannot 
        be used if you enable midpoint correction.
    simulate_pval : bool 
        Indicate whether to compute p-values by Monte Carlo simulation,
         in larger than 2 x 2 tables.
    replicate : int
        An integer specifying the number of replicates used in the MonteCarlo
        test.
    workspace : int
        An integer specifying the workspace size. Default value is 300. 
    attempt : int
        Number of attempts to try, if the workspace size is not enough. 
        On each attempt, the workspace size is doubled. 
    seed : int
        Random number to use as seed. If a seed isn't provided. 4 bytes will be
        read from os.urandom. If this fail, getrandbits of the random module 
        (with 32 random bits) will be used. In the particular case where both
        failed, the current time will be used

    Returns
    -------
    p_value : float
        The probability of a more extreme table, where 'extreme' is in a 
        probabilistic sense.

    Notes
    -----
    The calculated odds ratio is different from the one R uses. This scipy
    implementation returns the (more common) "unconditional Maximum
    Likelihood Estimate", while R uses the "conditional Maximum Likelihood
    Estimate".
    For tables with large numbers, the (inexact) chi-square test implemented
    in the function `chi2_contingency` can also be used.
    Examples
    --------
    Say we spend a few days counting whales and sharks in the Atlantic and
    Indian oceans. In the Atlantic ocean we find 8 whales and 1 shark, in the
    Indian ocean 2 whales, 5 sharks and in the Pacific 12 whales and 2 sharks.
    Then our contingency table is::
                Atlantic  Indian    Pacific
        whales     8        2       12
        sharks     1        5       2
    We use this table to find the p-value:
    >>> from Fisher import fisher_exact
    >>> pvalue = stats.fisher_exact([[8, 2, 12], [1, 5, 2]])
    >>> pvalue
    0.01183...
    """

    workspace = 2 * int(workspace / 2)
    # int32 is not enough for the algorithm
    c = np.asarray(table, dtype=np.int64)

    if len(c.shape) > 2:
        raise ValueError(
            "The input `table` should not have more than 2 dimension.")

    if np.any(np.asarray(c.shape) < 2):
        raise ValueError("The input `table` must be at least of shape (2, 2).")

    # We expect all values to be non-negative
    if np.any(c < 0):
        raise ValueError("All values in `table` must be nonnegative.")

    nr, nc = c.shape

    if (nr == 2 and nc == 2):
        # I'm not sure what the fisher_exact module of ss do.
        # So use my own function if midp is asked
        if not midP:
            # in this case, just use the default scipy
            # could remove this in the future
            return ss.fisher_exact(c, alternative)[1]
        else:
            return _midp(c)

    else:
        pval = None
        if simulate_pval:
            sr = c.sum(axis=1)
            sc = c.sum(axis=0)
            # The zero colums and rows are droped here, see R function
            c = c[sr > 0, :][:, sc > 0]
            nr, nc = c.shape
            if nr < 2 or nc < 2:
                raise ValueError(
                    'Less than 2 non-zero column or row marginal,\n %s' % c)

            statistic = -np.sum(lgamma(c + 1))
            tmp_res = _fisher_sim(c, replicate, seed)
            almost = 1 + 64 * np.finfo(np.double).eps
            pval = (1 + np.sum(tmp_res <= statistic / almost)) / \
                (replicate + 1.)
        elif hybrid:
            expect, percnt, emin = 5, 80, 1  # this is the cochran condition
            pval = _execute_fexact(nr, nc, c, nr, expect,
                                   percnt, emin, workspace, attempt, midP)
        else:
            expect, percnt, emin = -1, 100, 0
            pval = _execute_fexact(nr, nc, c, nr, expect,
                                   percnt, emin, workspace, attempt, midP)

        return pval


def _execute_fexact(nr, nc, c, nnr, expect, percnt, emin, workspace,
                    attempt=2, midP=False):
    """Execute fexact using the fortran routine"""

    pval = None
    success = False
    ntry = 0
    while not (success or ntry >= attempt):
        ntry += 1
        try:
            pval = f.fexact(nr, nc, c, nnr, expect, percnt, emin, workspace)
            success = True
        except Exception as e:
            logging.warning(
                "Workspace : %d is not enough. You should increase it.")
    if midP:
        return pval[1] - pval[0] * 0.5
    else:
        return pval[1]


def _fisher_sim(c, replicate, seed=None):
    """Performs a simulation with `replicate` replicates in order to find an 
    alternative contingency test with the same margin.
    Parameters
    ----------
    c : array_like of ints
        A m x n contingency table.  Elements should be non-negative integers.
    replicate : int
        Number of replicates to perform for the simulation

    seed : int
        A random number to be used as seed
    """
    if seed is None:
        try:
            seed = os.urandom(4)
            seed = int(seed.encode('hex'), 16)
        except:
            try:
                seed = int(random.getrandbits(32))
            except:
                import time
                seed = int(time.time())

    seed = np.array([seed], dtype='int32')
    key = np.array([False], dtype=bool)
    ierror = np.array([0], dtype='int32')
    sr, sc = c.sum(axis=1).astype('int32'), c.sum(axis=0).astype('int32')
    nr, nc = len(sr), len(sc)
    n = np.sum(sr)
    results = np.zeros(replicate)

    fact = np.zeros(n + 1)
    for i in xrange(2, n + 1):
        fact[i] = fact[i - 1] + np.log(i)

    observed = np.zeros((nr, nc), dtype="int32", order='F')
    for it in xrange(replicate):
        rcont2(nrow=nr, ncol=nc, nrowt=sr, ncolt=sc, key=key,
               seed=seed, matrix=observed, ierror=ierror)
        # if we do not have an error, make spcial action
        ans = 0.
        tmp_observed = observed.ravel()
        if ierror[0] != 0:
            raise ValueError("Fortran subroutine rcont2 return an error !")
        for j in xrange(nc):
            i = 0
            ii = j * nr
            while(i < nr):
                ans -= fact[tmp_observed[ii]]
                i += 1
                ii += 1
        results[it] = ans
    return results


def _midp(c):
    """Performs Fisher's Exact test with midp correction
    Parameters
    ----------
    c : array_like of ints
        A m x n contingency table. Elements should be non-negative integers. 
    """
    sr, sc = c.sum(axis=1).astype('int32'), c.sum(axis=0).astype('int32')
    nr, nc = len(sr), len(sc)
    n = np.sum(sr)
    global result
    result = []
    logfact = np.zeros(n + 1)
    for i in xrange(2, n + 1):
        logfact[i] = logfact[i - 1] + np.log(i)

    def callback(iflag, table, m, n, rowsum, colsum, prob, mult):
        global result
        if iflag == 2 or (isinstance(iflag, list) and iflag[0] == 2):
            tmp = np.asarray(table.tolist())
            if not np.all(c == tmp):
                # only add table that are different from the original that we
                # had
                result.append(tmp.ravel())

    def table_prob(table, sr, sc, n):
        num = np.sum([logfact[i] for i in sr]) + \
            np.sum([logfact[j] for j in sc])
        den = logfact[n] + np.sum([logfact[i] for i in table])
        return np.exp(num - den)

    # After this call, all the result should be in result
    rcont(sr, sc, callback, 0)
    problist = []
    for rtable in result:
        # compute the probability of each table
        problist.append(table_prob(rtable, sr, sc, n))

    problist = np.asarray(problist)
    problist.sort()
    cur_prob = table_prob(c.ravel(), sr, sc, n)
    pval = np.sum(problist[problist <= cur_prob])
    return pval + cur_prob * 0.5
