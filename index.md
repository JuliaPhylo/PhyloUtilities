---
layout: default
title: About
nav_order: 1
---

This site has material for a software workshop on
phylogenetic networks,
used at the 2018 & 2019 [MBL](https://molevol.mbl.edu/index.php/Main_Page)
workshop on molecular evolution (earlier version for a [2016 workshop](http://tandy.cs.illinois.edu/symposium-2016.html)).
<!-- also at the 2019 [Midwest](https://phylosdd.github.io/MidwestPhylo2019/) Phylogenetics Workshop -->
It covers steps of the TICR pipeline to go from a bunch of multiple alignments
(aligned gene sequences, or loci) to a concordance factors. These concordance factors can then be used to estimate
 a phylogenetic network that displays the relationships between the species in the alignments.

There is another newer [online tutorial](https://solislemuslab.github.io/snaq-tutorial/) created as part of the workshop in the Kew Royal Botanical Gardens: [Methodological Advances in Reticulate Evolution](https://gtiley.github.io/RBG-Networks/about/) taught in November 2023.

## topics covered

<!-- - [requirements](https://github.com/JuliaPhylo/v0.16PhyloNetworks-wiki-tutorial/wiki/Workshop-Requirements) -->
- [example data](https://github.com/JuliaPhylo/v0.16PhyloNetworks-wiki-tutorial/wiki/Example-Data)
  to download
- [TICR pipeline](https://github.com/JuliaPhylo/v0.16PhyloNetworks-wiki-tutorial/wiki/TICR-from-alignments-to-quartet-concordance-factors)
  overview:
  from sequences to quartet concordance factors
  (CFs, proportion of genes having a particular history)
- [TICR test](https://github.com/JuliaPhylo/v0.16PhyloNetworks-wiki-tutorial/wiki/TICR-test-tree-versus-network):
  is a population tree with ILS sufficient (vs network)?

## Set-up

### locally

- Download [BUCKy](http://pages.stat.wisc.edu/~ane/bucky/index.html)
- Download [TICR](https://github.com/nstenz/TICR)
- Download [QuartetMaxCut](http://research.haifa.ac.il/%7Essagi/software/QMCN.tar.gz)
- Download [MrBayes](http://nbisweden.github.io/MrBayes/)

### for MBL workshop

Login to your particular virtual machine (VM) using the IP address on the sticker attached to the back of your name tag. If, for example, your IP address was 123.456.789.321, you would type the following into the terminal on your local computer (i.e. your laptop) and press the enter key:

```bash
ssh moleuser@123.456.789.321
```

After login, you want to copy the `phylo-networks` in your home directory:

```
cp -r moledata/phylo-networks ./
cd phylo-networks
```

## more details

- [TICR](https://github.com/nstenz/TICR) pipeline

    - analyze each locus with MrBayes
    - do a concordance analysis with BUCKy on each set of 4 taxa
    - summarize all quartet concordance factors (CFs)
    - estimate a species tree using Quartet MaxCut

    We will *not* cover an alternative pipeline
    (which you could use outside the workshop) to

    - analyze each locus with RAxML, including bootstrap
    - estimate a species tree with ASTRAL


- TICR test in the R package [phylolm](https://github.com/lamho86/phylolm)
  (Testing Incongruence Checking in R)

   - test if a population tree with the coalescent (ILS only) is adequate
     to explain the quartet concordance factors
   - identify taxa involved in outlier quartets
