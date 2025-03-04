---
layout: default
title: About
nav_order: 1
---

# From alignments to concordance factors

This site covers the steps of the TICR pipeline to go from a bunch of multiple alignments
(aligned gene sequences, or loci) to a concordance factors. These concordance factors can then be used to estimate a phylogenetic network that displays the relationships between the species in the alignments.

- The original [TICR](https://github.com/nstenz/TICR) pipeline

    - analyze each locus with MrBayes
    - do a concordance analysis with BUCKy on each set of 4 taxa
    - summarize all quartet concordance factors (CFs)
    - estimate a species tree using Quartet MaxCut

    Or an alternative pipeline to

    - analyze each locus with RAxML, including bootstrap
    - estimate a species tree with ASTRAL


- TICR test in the R package [phylolm](https://github.com/lamho86/phylolm)
  (Testing Incongruence Checking in R), also now in Julia in [QuartetNetworkGoodnessFit.jl](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl)

   - test if a population tree with the coalescent (ILS only) is adequate
     to explain the quartet concordance factors
   - identify taxa involved in outlier quartets

## Topics covered

- [Example data](https://juliaphylo.github.io/PhyloUtilities/notebooks/Example-Data.html)
- [Gene trees with MrBayes](https://juliaphylo.github.io/PhyloUtilities/notebooks/Gene-Trees-MrBayes.html)
- [Quartet CFs with BUCKy](https://juliaphylo.github.io/PhyloUtilities/notebooks/Quartet-CF-BUCKy.html)
- [Species tree with TreeQMC](https://juliaphylo.github.io/PhyloUtilities/notebooks/Species-tree-from-quartet-CFs-QMC.html)
- [Alternative pipeline: RAxML+ASTRAL](https://juliaphylo.github.io/PhyloUtilities/notebooks/Gene-Trees-RAxML.html)
- [TICR goodness-of-fit test](https://juliaphylo.github.io/PhyloUtilities/notebooks/TICR-test-tree-versus-network.html): is a population tree with ILS sufficient (vs network)?

{: .note }
**Note** that all these sections are covering the original TICR pipeline well suited for a machine or a cluster of machines without a job scheduler. The scripts automatically parallelize the work across the available cores. Users with access to job schedulers like [SLURM](https://slurm.schedmd.com/) or
[SGE](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) should check the section on [SLURM pipeline](). The job scheduler does the work of parallelizing the work across available cores.
The scripts, in this second pipeline, were created to take full advantage
of job scheduler capabilities. They were developed for a cluster running SLURM. Adjustments to the submit scripts will be needed, to adapt to your own
SLURM configuration or to the syntax that your job scheduler wants.

## Set-up

### Locally

- Download [BUCKy](http://pages.stat.wisc.edu/~ane/bucky/index.html)
- Download [TreeQMC](https://github.com/molloy-lab/TREE-QMC)
- Download [MrBayes](http://nbisweden.github.io/MrBayes/)
- Git clone this repository: `git clone https://github.com/JuliaPhylo/PhyloUtilities.git`

### For MBL workshop

Login to your particular virtual machine (VM) using the IP address on the sticker attached to the back of your name tag. If, for example, your IP address was 123.456.789.321, you would type the following into the terminal on your local computer (i.e. your laptop) and press the enter key:

```bash
ssh moleuser@123.456.789.321
```

After login, you want to copy the `phylo-networks` in your home directory:

```
cp -r moledata/phylo-networks ./
cd phylo-networks
```

We will *not* cover an alternative RAxML+ASTRAL pipeline (which you could cover outside the workshop).


