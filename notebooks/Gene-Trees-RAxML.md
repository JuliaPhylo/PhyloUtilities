---
layout: default
title: RAxML+ASTRAL
nav_order: 6
---

# Estimating gene trees with RAxML and species tree with ASTRAL

A faster alternative to running MrBayes on each gene and BUCKy is to run instead
RAxML on each gene (with bootstrap for later use), followed by counting the number
of gene trees that have each quartet (which will be done within PhyloNetworks).

**For MBL workshop, this section is optional.**

The advantage of MrBayes+BUCKy is that uncertainty in gene trees is reduced and
integrated out to estimate quartet concordance factors. But it's slower.

To run RAxML on each gene, we can use the script `raxml.pl` in the `script/` folder.
To do so, still from the data folder `baseline.gamma0.3_n30/`, we first uncompress
the tarball that has the alignments in nexus format. The script won't do this for us
unfortunately. Also, to keep (thousands of) files organized, we uncompress these files
in a folder `nexus/` within our `input/` folder. Note that if you followed the steps in [The Data](https://github.com/JuliaPhylo/v0.16PhyloNetworks-wiki-tutorial/wiki/The-Data) section, you already created the `input/nexus` folder.

```
mkdir input/nexus
tar -xzf input/1_seqgen.tar.gz -C input/nexus
```

We can then run the script, which is not parallelized by the way and cannot send jobs to a
cluster (unlike the scripts from TICR).
RAxML sends a lot of output to the screen, 30 times if we have 30 genes, so we
redirect this overwhelming output (and any possible error message) to a file `mylog` like this:
```
../../scripts/raxml.pl --seqdir=input/nexus --raxmldir=raxml --astraldir=astral > mylog 2>&1 &
```
Let's check that everything ran smoothly and that we have our desired output.
We have a new `raxml/` directory with a bunch of things, including one file containing
the best tree from each gene:
```
$ ls raxml/
besttrees.tgz	besttrees.tre	bootstrap	contrees.tgz	raxml.pl.log

$ head -n 4 raxml/besttrees.tre
(((4:0.01015732493562751526,3:0.01431694705110660507):0.01028775678308786225,(1:0.00197193126842715389,2:0.00810229112284156609):0.01021762708823881899):0.02482468034550755487,5:0.03813009647365416671,6:0.07829158619530514340):0.0;
((3:0.01822494774942512788,(5:0.02442835505708105398,4:0.01822832410801497258):0.00000100000050002909):0.03638839679233905194,(2:0.00571736034091331266,1:0.00451901895337002840):0.04118871647981920542,6:0.08034263217082675268):0.0;
(((4:0.00606576976462778941,3:0.00607447630066324223):0.00338289128901973421,5:0.01486951087863153456):0.04130627869590708379,(1:0.00000100000050002909,2:0.00601218297756087993):0.04083936720011970695,6:0.05249779132102743578):0.0;
(6:0.05336929993901522173,((3:0.00402917441342817879,4:0.00403066318961684770):0.00659185294481937686,(1:0.00403034039908813698,2:0.00203303117370299017):0.01635211316824565497):0.02577384948609293819,5:0.05669026111440958471):0.0;
```
We also have a new folder `raxml/boostrap/` containing one bootstrap tree file per gene:
```
$ ls raxml/bootstrap/ | head
RAxML_bootstrap.1_seqgen1
RAxML_bootstrap.1_seqgen10
RAxML_bootstrap.1_seqgen11
RAxML_bootstrap.1_seqgen12
RAxML_bootstrap.1_seqgen13
RAxML_bootstrap.1_seqgen14
RAxML_bootstrap.1_seqgen15
RAxML_bootstrap.1_seqgen16
RAxML_bootstrap.1_seqgen17
RAxML_bootstrap.1_seqgen18
```
and we have a new folder `astral/` that has the results from ASTRAL (which the script runs
after RAxML by default):
```
$ ls astral/
BSlistfiles		astral.screenlog	astral.tre
```

Note: the `astral` command in the script
**assumes that there is a single individual per species**.

If your data has multiple alleles or individual per species, ASTRAL will need
to be given a mapping file, to know which individuals map to which species.
In this case, we can run `raxml.pl` without the final `astral` step (option `--nodoastral`), 
and then run ASTRAL separatedly with a specific mapping file.

That is, if your data has multiple alleles, you would run `raxml.pl` without the 
final `astral` step:
```shell
../../scripts/raxml.pl --seqdir=input/nexus --raxmldir=raxml --nodoastral > mylog 2>&1 &
```

The `raxml.pl` script will print at the end of its log file the specific ASTRAL command
(which has the appropriate path to the input files for ASTRAL) that you would need to run next:
```bash
$ tail -n 2 raxml/raxml.pl.log
astral could be run with:
java -jar astral -i raxml/besttrees.tre -b astral/BSlistfiles -r 100 -o astral/astral.tre > astral/astral.screenlog 2>&1
```
You will need to modify this command by adding to it the option that gives the 
mapping file.

In our example, we had a single individual per species, so the ASTRAL command
that the script ran was just what we needed. In the output `.tre` file,
we have 100 bootstrap species trees, and another 2 extra trees: the original
species tree (on the original gene trees) twice: first with edges annotated
with bootstrap values, second with edges annotated with bootstrap values and
edge lengths (in coalescent units).

```shell
$ head -3 astral/astral.tre
(3,(4,((5,6)1:0.5764721782579237,(2,1)1:2.483246898369636)1:1.154550029033082)); 
(3,(4,((5,6)1:0.6311555045319776,(2,1)1:2.68060633252813)1:1.1234191104379085)); 
(3,(4,((5,6)1:0.6172017897581119,(2,1)1:2.68060633252813)1:1.1157854855828375)); 

$ tail -n 5 astral/astral.tre 
(3,(4,((5,6)1:0.6405673387143238,(2,1)1:2.610402073854883)1:1.1082096917743796)); 
(3,(4,((5,6)1:0.5545414747639535,(2,1)1:2.68060633252813)1:1.0932281381587627)); 
(3,(4,((5,6)1:0.5898652355943617,(2,1)1:2.7561138850362754)1:1.1006908593603526)); 
(3,(4,((5,6)100.0,(1,2)100.0)100.0)); 
(3,(4,((5,6)100.0:0.6548532959618005,(2,1)100.0:2.6448882499260518)100.0:1.1388634328653822)); 
```
In this `astral/` folder, we also have the list of all bootstrap RAxML files, which was needed
for ASTRAL, and which will be reused in PhyloNetworks to get bootstrap species networks:
```shell
$ head astral/BSlistfiles 
raxml/bootstrap/RAxML_bootstrap.1_seqgen1
raxml/bootstrap/RAxML_bootstrap.1_seqgen10
raxml/bootstrap/RAxML_bootstrap.1_seqgen11
raxml/bootstrap/RAxML_bootstrap.1_seqgen12
raxml/bootstrap/RAxML_bootstrap.1_seqgen13
raxml/bootstrap/RAxML_bootstrap.1_seqgen14
raxml/bootstrap/RAxML_bootstrap.1_seqgen15
raxml/bootstrap/RAxML_bootstrap.1_seqgen16
raxml/bootstrap/RAxML_bootstrap.1_seqgen17
raxml/bootstrap/RAxML_bootstrap.1_seqgen18
```
Now that we got all our results with no error, we can remove the uninteresting
output from RAxML in `mylog`:

    rm mylog



