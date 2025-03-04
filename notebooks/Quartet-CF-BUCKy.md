---
layout: default
title: Quartet CFs with BUCKy
nav_order: 4
---

# Estimating quartet CFs with BUCKy

We can combine the MrBayes analyses across all loci to estimate the proportion of genes
(concordance factor) having a given quartet, for all quartets:
```
$ ../../scripts/bucky.pl mb-output/1_seqgen.mb.tar -o bucky-output

Checking for BUCKy version >= 1.4.4...
  BUCKy version: 1.4.4.
  BUCKy version check passed.

Script was called as follows:
perl bucky.pl mb-output/1_seqgen.mb.tar -o bucky-output

Found 6 taxa shared across all genes in this archive, 15 of 15 possible quartets will be run using output from 30 total genes.
Summarizing MrBayes output for 30 genes.
Job server successfully created.

  Analyses complete: 15/15.
  All connections closed.
Total execution time: 57 seconds.
```
The script created a new directory named `bucky-output` (like we asked above), which contains
several tarballs (results from MrBayes, mbsum and from bucky)
```
$ ls
astral	bucky-output  input  mb-output	raxml  snaq

$ ls bucky-output/
1_seqgen.BUCKy.tar  1_seqgen.mb.tar
1_seqgen.CFs.csv    1_seqgen.mbsum.tar.gz
```
The `.csv` file (spreadsheet) lists all the 4-taxon sets and their estimated quartet
concordance factors:
```
$ column -s, -t bucky-output/1_seqgen.CFs.csv
taxon1  taxon2  taxon3  taxon4  CF12.34               CF12.34_lo          CF12.34_hi          CF13.24               CF13.24_lo          CF13.24_hi          CF14.23               CF14.23_lo         CF14.23_hi         ngenes
1       3       5       6       0.565033333333333     0.5                 0.633333333333333   0.0903                0.0666666666666667  0.133333333333333   0.344666666666667     0.3                0.4                30
1       3       5       4       0.0005                0                   0                   0.8599                0.833333333333333   0.866666666666667   0.139566666666667     0.133333333333333  0.166666666666667  30
2       3       5       4       0.0005                0                   0                   0.8599                0.833333333333333   0.866666666666667   0.139566666666667     0.133333333333333  0.166666666666667  30
2       3       5       6       0.565033333333333     0.5                 0.633333333333333   0.0903                0.0666666666666667  0.133333333333333   0.344666666666667     0.3                0.4                30
1       3       6       4       3.33333333333333e-05  0                   0                   0.8885                0.866666666666667   0.9                 0.111466666666667     0.1                0.133333333333333  30
2       1       5       6       0.999866666666667     1                   1                   6.66666666666667e-05  0                   0                   6.66666666666667e-05  0                  0                  30
2       5       6       4       0.0401666666666667    0.0333333333333333  0.0666666666666667  0.263066666666667     0.233333333333333   0.333333333333333   0.696766666666667     0.633333333333333  0.733333333333333  30
2       1       3       5       1                     1                   1                   0                     0                   0                   0                     0                  0                  30
2       1       6       4       0.999866666666667     1                   1                   6.66666666666667e-05  0                   0                   6.66666666666667e-05  0                  0                  30
2       1       5       4       1                     1                   1                   0                     0                   0                   0                     0                  0                  30
1       5       6       4       0.0401666666666667    0.0333333333333333  0.0666666666666667  0.263066666666667     0.233333333333333   0.333333333333333   0.696766666666667     0.633333333333333  0.733333333333333  30
2       1       3       6       0.999866666666667     1                   1                   6.66666666666667e-05  0                   0                   6.66666666666667e-05  0                  0                  30
2       1       3       4       1                     1                   1                   0                     0                   0                   0                     0                  0                  30
3       5       6       4       0.0731666666666667    0.0666666666666667  0.1                 0.0424666666666667    0.0333333333333333  0.0666666666666667  0.884366666666667     0.833333333333333  0.9                30
2       3       6       4       3.33333333333333e-05  0                   0                   0.8885                0.866666666666667   0.9                 0.111466666666667     0.1                0.133333333333333  30
```
<!-- for 300 genes: <br><img src="screenshots/CFtable.png" width="800"> -->

This will be the input to
[SNaQ](https://juliaphylo.github.io/SNaQ.jl/stable/):
to estimate a phylogenetic tree or network.
SNaQ will also need a topology to start the network search,
like a first estimate of a species tree. We can use a tree from TreeQMC (next step).
Alternatively, we can run RAxml on each locus and then run ASTRAL
[(here)](https://juliaphylo.github.io/PhyloUtilities/notebooks/Gene-Trees-RAxML.html).

Back to BUCKy analyses: If the tree samples from MrBayes have already been summarized with
`mbsum`, the script `bucky.pl` can take these in (instead of tree samples from MrBayes)
using the option `-s`. In this case, the output files from `mbsum` need to be bundled together,
like this for instance: `tar czf my-mbsum.tar.gz *.in` if the `mbsum` files end with `.in`.
An example of such bundle is `bucky-output/1_seqgen.mbsum.tar.gz`.
Then run `bucky.pl` like this:
`bucky.pl my-mbsum.tar.gz -s -o bucky-output`.
