---
layout: default
title: Example Data
nav_order: 2
---

There are two ways to download the input data files:

- navigate to where you want the data, then clone the repository that
  contains both this tutorial and the files:

  `git clone https://github.com/JuliaPhylo/v0.16PhyloNetworks-wiki-tutorial.wiki.git`

  The repository will appear as a directory named `PhyloNetworks.jl.wiki`.
  Inside, you should find a directory `data_results`:

  ```
  cd PhyloNetworks.jl.wiki
  ls
  cd data_results
  ls
  ```

- or download and uncompress a tarball:
  [.zip](http://www.stat.wisc.edu/~ane/PhyloNetworks/data_results.zip)
  or [.tgz](http://www.stat.wisc.edu/~ane/PhyloNetworks/data_results.tgz)

These tutorial files are organized in a directory `data_results/` that has:
- 6 example data sets, with a separate directory for each.
  These are all simulated data, for which we know the true network.

  - 1 data set has 15 taxa (`n15.gamma0.30.20.2_n300`).
  - 5 data sets have 6 taxa (`baseline.*` and `n6.gamma0.30.2_n300`), 
    to run fast during the tutorial. Among those,
    4 (`baseline.*`) were generated under a network with 1 reticulation and different numbers of genes
    (from 30 to 1000 genes), and 1 data set was generated under a network with 2 reticulations
    (`n6.gamma0.30.2_n300`).

  These examples will provide ways to look at the effect of the number of taxa, the number
  of genes, and the number of reticulations.

- a folder `scripts/` containing all the scripts we will use in the tutorial.


Within each dataset folder, you will find several directories:

- `input/`: contains input data. One file has the true gene trees (simulated from a network)
  and the other file is a tarball containing all the alignments (e.g. 1000 alignments if the
  data set has 1000 genes).
  All the other folders and files can be recreated using the scripts. The main results files
  are provided, though, to allow participants to pick up the tutorial at any step.
- `bucky-output/`: contains the main output file from running MrBayes + BUCKy, which is a table
  listing quartet concordance factors (CFs). Also contains a file with the species estimated
  from quartet CFs using Quartet MaxCut.
- `raxml/`: contains one file with all the best RAxML trees, one for each gene;
  and a directory `bootstrap/` containing one file with 100 bootstrap for each gene.
- `astral/`: contains one file that lists all the bootstrap tree files from RAxML, and
  another file with the results of ASTRAL. These results consists of 102 trees: first 100
  bootstrap species trees, then their consensus, then the species tree estimated from the
  original data annotated with bootstrap support.
- `snaq/`: contains various files for the various estimated networks: best networks
  and bootstrap networks.

For 4 of the 6 data sets, the main results of each step are provided. Participants can
continue on the next step even if the previous step did not work on their laptop.
