Module Fisher
-------------

Variables
---------
f

lgamma

rcont2

Functions
---------
execute_fexact(nr, nc, c, nnr, expect, percnt, emin, workspace, attempt=2)
    Execute fexact using the fortran routine

fisher_exact(table, alternative='two-sided', hybrid=False, simulate_pval=False, replicate=2000, workspace=300, attempt=2, seed=None)
    Performs a Fisher exact test on a 2x2 contingency table.
    Parameters
    ----------
    table : array_like of ints
        A 2x2 contingency table.  Elements should be non-negative integers.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Which alternative hypothesis to the null hypothesis the test uses.
        Default is 'two-sided'.  Only used in the 2 x 2 case.
    mult : int 
        Specify the size of the workspace used in the network algorithm.  
        Only used for non-simulated p-values larger than 2 x 2 table. 
        You might want to increase this if the p-value failed
    hybrid : bool
        Only used for larger than 2 x 2 tables, in which cases it indicates
        whether the exact probabilities (default) or a hybrid approximation 
        thereof should be computed.
    simulate_pval : bool 
        Indicate whether to compute p-values by Monte Carlo simulation,
         in larger than 2 x 2 tables.
    replicate : int
        An integer specifying the number of replicates used in the Monte Carlo test.
    workspace : int
        An integer specifying the workspace size. Default value is 300. 
    attempt : int
        Number of attempts to try, if the workspace size is not enough. 
        On each attempt, the workspace size is doubled. 
    seed : int
        Random number to use as seed. If a seed isn't provided. 4 bytes will be read 
        from os.urandom. If this fail, getrandbits of the random module 
        (with 32 random bits) will be used. In the particular case where both failed,
        the current time will be used

    Returns
    -------
    prt : float
        Probability of the observed table for fixed marginal totals.
        If it's a 2x2 table, the prior odds ratio is returned instead.
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
    Indian ocean 2 whales and 5 sharks. Then our contingency table is::
                Atlantic  Indian    Pacific
        whales     8        2       12
        sharks     1        5       2
    We use this table to find the p-value:
    >>> from Fisher import fisher_exact
    >>> oddsratio, pvalue = stats.fisher_exact([[8, 2, 12], [1, 5, 2]])
    >>> pvalue
    0.01183...

fisher_sim(c, replicate, seed=None)
    Performs a simulation with `replicate` replicates in order to find an 
     alternative contingency test with the same margin.
    Parameters
    ----------
    c : array_like of ints
        A m x n contingency table.  Elements should be non-negative integers.
    replicate : int
        Number of replicates to perform for the simulation
        
    seed : int
        A random number to be used as seed
