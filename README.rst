FisherExact
===========

Fisher exact test for mxn contingency table

Installation
------------

You can install it with pip : pip install FisherExact

Binary Usage
------------

A binary is provided to use FisherExact from the terminal

::

    usage: fexact [-h] [--simulate [SIMULATE]] [--hybrid] [--midP]
                  [--retry ATTEMPT] [--workspace WORKSPACE] [--version]
                  table

    Fisher's Exact test for mxn contingency table

    positional arguments:
      table                 Contingency table in a file, without header

    optional arguments:
      -h, --help            show this help message and exit
      --simulate [SIMULATE]
                            Simulate p-values with n replicates
      --hybrid              Use hybrid mode
      --midP                Use midP correction
      --retry ATTEMPT       Number of attempt to made if execution fail
      --workspace WORKSPACE
                            Workspace size to use, Increase this if the program
                            crash
      --version             show program's version number and exit

Contingency table format if fexact is used as binary
----------------------------------------------------

The accepted format is space/tab or comma (or both) separated values
with an optionnal first line starting with a "#" that specified the
number of rows and column:

For example, the following format are accepted

::

    # 2 3
    8        2       12
    1        5       2

::

    8 2 12
    1 5 2

::

    #2, 3
    8   2   12
    1   5   2

::

    8,2,12
    1,5,2


Functions
---------
    fisher_exact(table, alternative='two-sided', hybrid=False, midP=False, simulate_pval=False, replicate=2000, workspace=300, attempt=2, seed=None)
        Performs a Fisher exact test on a 2x2 contingency table.
        Parameters
        ----------
        table : array_like of ints
            A mxn contingency table.  Elements should be non-negative integers.
        alternative : {'two-sided', 'less', 'greater'}, optional
            Which alternative hypothesis to the null hypothesis the test uses.
            Default is 'two-sided'.  Only used in the 2 x 2 case (with the scipy function).
            In every other case, the two-sided pval is returned.
        mult : int 
            Specify the size of the workspace used in the network algorithm.  
            Only used for non-simulated p-values larger than 2 x 2 table. 
            You might want to increase this if the p-value failed
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
        Indian ocean 2 whales and 5 sharks and in the Pacific 12 whales andd 2 sharks
        Then our contingency table is::
                    Atlantic  Indian    Pacific
            whales     8        2       12
            sharks     1        5       2
        We use this table to find the p-value:
        >>> from Fisher import fisher_exact
        >>> pvalue = stats.fisher_exact([[8, 2, 12], [1, 5, 2]])
        >>> pvalue
        0.01183...