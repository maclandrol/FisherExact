# FisherExact

Fisher exact test for mxn contingency table

## Installation 

You can install it with pip : pip install FisherExact

## Binary Usage

A binary is provided to use FisherExact from the terminal

	```usage: fexact [-h] [--simulate [SIMULATE]] [--hybrid] [--midP]
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
	  --version             show program's version number and exit```


## Contingency table format if fexact is used as binary

The accepted format is space/tab or comma (or both) separated values with an optionnal first line starting with a "#" that specified the number of rows and column:

For example, the following format are accepted

```
# 2 3
8        2       12
1        5       2
```

```
8 2 12
1 5 2
```

```
#2, 3
8	2	12
1	5	2
```

```
8,2,12
1,5,2
```


