---
layout: default
title: Species tree from quartet CFs
nav_order: 5
---

# Estimating a species tree from quartet CFs with QMC

Quartet MaxCut [(QMC)](https://pubmed.ncbi.nlm.nih.gov/21030737/)
uses quartets as input, and returns a tree containing as many input quartets as possible.
The goal is similar to the goal of ASTRAL, but the input is different:
quartets instead of gene trees. We will use this QMC tree as a starting tree
to search for a better tree (or network) under a likelihood-based criterion.

To extract quartets with highest CFs and use them as input to QMC to get a tree,
we can use the `get-pop-tree.pl` script:

```
$ ../../scripts/get-pop-tree.pl bucky-output/1_seqgen.CFs.csv

Script was called as follows:
perl get-pop-tree.pl bucky-output/1_seqgen.CFs.csv

Parsing major resolution of each 4-taxon set... done.
Running Quartet Max Cut...

Quartet MaxCut version 2.10 by Sagi Snir, University of Haifa
quartet file is 1_seqgen.QMC.txt,
Number of quartets is 15, max element 6
Number of quartets read: 15, max ele 6
Started working at  Mon May 23 09:44:12 2016
Ended working at  Mon May 23 09:44:12 2016

Quartet Max Cut complete, tree located in '1_seqgen.QMC.tre'.
```

We have a new file `1_seqgen.QMC.tre` with our QMC tree.
And can look at it, then organize results.
A good place this tree in the bucky folder:

    $ cat 1_seqgen.QMC.tre
    $ mv 1_seqgen.QMC.tre bucky-output
    $ ls bucky-output/
    1_seqgen.BUCKy.tar 1_seqgen.CFs.csv 1_seqgen.QMC.tre 1_seqgen.mb.tar 1_seqgen.mbsum.tar.gz


