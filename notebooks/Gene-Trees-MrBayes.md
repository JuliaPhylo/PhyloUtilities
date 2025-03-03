---
layout: default
title: Gene trees with MrBayes
nav_order: 3
---

# Estimating gene trees with MrBayes

We want to analyze each of the 30 loci with MrBayes.
First, make sure you have MrBayes installed.
On the MBL cluster:

```
$ which mb
/usr/local/bin/mb
```

Next, choose settings for MrBayes: model, prior for branch lengths etc.
Save them in a MrBayes block. Below: HKY model, 100,000 generations
2 chains (1 cold & 1 heated), 2 independent runs.
These settings were chosen to run things fast during this tutorial, but for a
real data set different setting should be chosen (such as 1 million generations,
3 chains and 3 runs).

```
$ cat ../../scripts/mbblock.txt
begin mrbayes;
set nowarnings=yes;
set autoclose=yes;
lset nst=2;
mcmcp ngen=100000 burninfrac=.25 samplefreq=50 printfreq=10000 [increase these for real]
diagnfreq=10000 nruns=2 nchains=2 temp=0.40 swapfreq=10;       [increase for real analysis]
mcmc;
sumt;
end;
```

We are ready to analyze all loci with MrBayes:
```
$ ../../scripts/mb.pl input/1_seqgen.tar.gz -m ../../scripts/mbblock.txt -o mb-output

Script was called as follows:
perl mb.pl input/1_seqgen.tar.gz -m ../../scripts/mbblock.txt -o mb-output

Appending MrBayes block to each gene... done.

Job server successfully created.

  Analyses complete: 30/30.
  All connections closed.
Total execution time: 46 seconds.
```

If a cluster is available with different machines, analyses can be parallelized
across machines (not just across nodes of the same machine) by adding an option
`--machine-file hosts.txt`, where `hosts.txt` is a simple text
file listing the machines available to use, in the format `user_name@machine_address`.
This file might look like this:
<br><img src="screenshots/hosts.png" width="300">


The script created a new directory named `mb-output` (like we asked above),
which contains a compressed tarball of all MrBayes output: `mb-output/1_seqgen.mb.tar`
```
$ ls
input  mb-output

$ ls mb-output/
1_seqgen.mb.tar	1_seqgen.tar.gz

$ tar -tf mb-output/1_seqgen.mb.tar
1_seqgen12.nex.tar.gz
1_seqgen11.nex.tar.gz
1_seqgen10.nex.tar.gz
...
1_seqgen7.nex.tar.gz
1_seqgen8.nex.tar.gz
1_seqgen9.nex.tar.gz
```

Let's look at the result file for the first locus. For this, let's go into
the new `mb-output` folder (we will need to go back to the main folder later),
create a folder `1_seqgen.mb` for exanding all the MrBayes results,
decompress the results for the first locus, and finally look at them:

```bash
cd mb-output
mkdir 1_seqgen.mb
tar -xvf 1_seqgen.mb.tar -C 1_seqgen.mb
ls 1_seqgen.mb
cd 1_seqgen.mb
mkdir 1_seqgen1.nex
tar -xzvf 1_seqgen1.nex.tar.gz -C 1_seqgen1.nex
ls 1_seqgen1.nex
```

We find a bunch of output including the log from MrBayes
(useful to track down bugs, if any) and the sample of
trees from each run (`*.t`), which will serve as input for BUCKy.
```
$ ls 1_seqgen1.nex
1_seqgen1.nex.ckp      1_seqgen1.nex.mcmc    1_seqgen1.nex.run2.t   1_seqgen1.nex.vstat
1_seqgen1.nex.ckp~     1_seqgen1.nex.parts   1_seqgen1.nex.run2.p
1_seqgen1.nex.con.tre  1_seqgen1.nex.run1.p  1_seqgen1.nex.trprobs
1_seqgen1.nex.log      1_seqgen1.nex.run1.t  1_seqgen1.nex.tstat
```

The most important files are those containing the sampled gene trees:

```
$ less -S 1_seqgen1.nex/1_seqgen1.nex.run1.t

#NEXUS
[ID: 5353870756]
[Param: tree]
begin trees;
   translate
       1 6,
       2 5,
       3 1,
       4 2,
       5 3,
       6 4;
   tree gen.0 = [&U] ((4:2.000000e-02,(6:2.000000e-02,2:2.000000e-02):2.000000e-02):2.000000e-02,(5:2.000000e-02,3:2.000000e-02):2.000000e-02,1:2.000000e-02);
   tree gen.50 = [&U] (((4:6.749905e-03,3:1.396825e-02):1.069614e-02,(5:2.164948e-02,6:7.205712e-03):1.861635e-02):3.204673e-02,2:2.699199e-02,1:3.668708e-02);
   tree gen.100 = [&U] (((4:5.960128e-03,3:6.373705e-03):1.170556e-02,(5:1.603174e-02,6:6.627496e-03):8.926058e-03):2.721952e-02,2:2.638575e-02,1:6.115106e-02);
   tree gen.150 = [&U] (((4:9.868770e-03,3:2.805685e-03):1.086845e-02,(5:2.218544e-02,6:6.324761e-03):9.347534e-03):1.152684e-02,2:5.165054e-02,1:6.528826e-02);
   tree gen.200 = [&U] (((4:1.098244e-02,3:2.853569e-03):1.351000e-02,(5:1.345820e-02,6:1.902973e-02):1.005224e-02):1.344231e-02,2:4.246484e-02,1:8.722215e-02);
   ...
   tree gen.99850 = [&U] (((4:2.183261e-02,3:4.841606e-03):4.529830e-03,(5:1.139391e-02,6:1.407128e-02):1.184515e-02):3.180385e-02,2:2.979942e-02,1:5.423444e-02)
   tree gen.99900 = [&U] (((6:1.200251e-02,5:1.781890e-02):1.025264e-02,(4:1.707598e-02,3:2.540760e-04):1.734718e-02):2.143296e-02,2:2.929320e-02,1:7.620077e-02)
   tree gen.99950 = [&U] (((4:2.988533e-03,3:3.844257e-03):1.021980e-02,(5:2.239675e-02,6:1.117266e-02):8.290780e-03):2.668343e-02,2:2.229218e-02,1:4.864512e-02)
   tree gen.100000 = [&U] (((6:1.238475e-02,5:1.200429e-02):7.875171e-03,(4:5.683504e-03,3:3.564505e-03):1.675113e-02):2.143296e-02,2:2.929320e-02,1:7.620077e-02
end;
```

Type `G` to go to the end of the file, `g` to come back to the beginning,
and `q` to quit viewing the file.

There are also trees from the second independent run:
```
less -S 1_seqgen1.nex/1_seqgen1.nex.run2.t
```

We won't be using it (because we will be using the full list of all sampled trees),
but we do have 50% majority rule consensus tree from the posterior distribution:

```
$ cat 1_seqgen1.nex/1_seqgen1.nex.con.tre
#NEXUS
[ID: 5353870756]
begin taxa;
	dimensions ntax=6;
	taxlabels
		6
		5
		1
		2
		3
		4
		;
end;
begin trees;
	translate
		1	6,
		2	5,
		3	1,
		4	2,
		5	3,
		6	4
		;
   tree con_50_majrule = [&U] (1[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:7.315541e-02[&length_mean=7.35856944e-02,length_median=7.31554100e-02,length_95%HPD={5.05359800e-02,9.70821400e-02}],2[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.569409e-02[&length_mean=3.66608383e-02,length_median=3.56940900e-02,length_95%HPD={2.11883100e-02,5.53439400e-02}],((3[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.064436e-03[&length_mean=3.66163042e-03,length_median=3.06443600e-03,length_95%HPD={8.10284600e-06,8.60207700e-03}],4[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:8.642618e-03[&length_mean=9.32276055e-03,length_median=8.64261800e-03,length_95%HPD={2.03317100e-03,1.75520400e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.098585e-02[&length_mean=1.16285593e-02,length_median=1.09858500e-02,length_95%HPD={3.41401000e-03,2.08909000e-02}],(5[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.455524e-02[&length_mean=1.51849277e-02,length_median=1.45552400e-02,length_95%HPD={5.76476500e-03,2.61160100e-02}],6[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.063229e-02[&length_mean=1.13693922e-02,length_median=1.06322900e-02,length_95%HPD={3.20096100e-03,2.05437800e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.101159e-02[&length_mean=1.16580154e-02,length_median=1.10115900e-02,length_95%HPD={3.05337700e-03,2.09244800e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:2.410702e-02[&length_mean=2.48581882e-02,length_median=2.41070200e-02,length_95%HPD={1.12175500e-02,3.92761300e-02}]);
end;
```

Before we move on to combine all these individual-locus analyses,
let's navigate back to the main folder for our data:
```
$ cd ../..
$ pwd
/home/moleuser/phylo-networks/data_results/baseline.gamma0.3_n30
```
