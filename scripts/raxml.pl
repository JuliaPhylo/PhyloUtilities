#!/usr/bin/perl

# May 26, 2016
# Perl script to run raxml + bootstrap for each locus.
# this is *not* part of the TICR pipeline
# adapted from Solis-Lemus, Yang and Ane (2016) scripts:
# https://github.com/crsl4/InconsistencySpeciesTreeGeneFlow/blob/master/scripts/estGeneTrees/raxml.pl

# model: HKY + Gamma (for rate variation across sites)

# usage:
#
# raxml.pl --seqdir=xxx/yyy --raxmldir=xxx --astraldir=xxx
#
# other options:
# --numCores = number of cores to use by RAxML. default 6
# --boot (default) or --noboot (not implemented!) to do RAxML bootstrap on each locus
# --numboot = number of bootstrap reps. default 100
# --conver2phylip (default) or --noconvert2phylip to convert nexus input files to phylip
# --doastral (default) or --nodoastral to do or not do ASTRAL at the end.
#
# the script will create a log file in raxmldir/raxml.pl.log

# warning: assuming all paths (seqdir and raxmldir) are relative paths, and linux/mac

use Getopt::Long;
use File::Path qw( make_path );
use strict;
use warnings;
use Carp;
#use lib '/u/c/l/claudia/lib/perl/lib/site_perl/'; # we need this because I had to install locally the Statistics module
#use Statistics::R;


# ================= parameters ======================

my $boot = 1; # boot=0 is not implemented, in fact!
my $numboot = 100;
my $numCores = 6;
my $seqdir;   # directory where sequences are
my $phylipdir;
my $raxmldir; # directory for output, including log for this script
my $astraldir;
my $convertphylip = 1;
my $doastral = 1;
my $raxml = 'raxmlHPC-PTHREADS-AVX'; # executable
my $astral = '/class/molevol-software/astral-5.6.3/astral.5.6.3.jar'; # adapt to your system

# -------------- read arguments from command-line -----------------------
GetOptions( 'numboot=i' => \$numboot,
	    'boot!' => \$boot,
	    'numCores=i' => \$numCores,
	    'seqdir=s' => \$seqdir,
	    'raxmldir=s' => \$raxmldir,
	    'astraldir=s' => \$astraldir,
	    'convert2phylip!' => \$convertphylip,
	    'doastral!' => \$doastral,
    );

die "seqdir not defined or not a directory" if (!(defined $seqdir) or !(-d $seqdir));
# sequence directory should have been uncompressed. If not, add this:
# system("tar -zxvf $seqdir");
if ($convertphylip) {
    $phylipdir = $seqdir;
    $phylipdir =~ s/[^\/]+$//; # remove last directory
    $phylipdir .= "phylip";    # replace by 'phylip'
    print "directory for phylip files: $phylipdir\n";
    make_path $phylipdir unless(-d $phylipdir);
}

die "a directory for RAxML output should be specified with the --raxmldir option" if (not defined $raxmldir);
die ("raxmldir should be only 1 level up\n") if ($raxmldir =~ /\//);
make_path $raxmldir unless(-d $raxmldir);
my $logfile = "$raxmldir/raxml.pl.log";

die "a directory for ASTRAL output should be specified with the --astraldir option" if ($doastral and (not defined $astraldir));
$astraldir = "astral" if !defined($astraldir);
make_path $astraldir if !(-d $astraldir);

my $currentdir = `pwd`;
chomp $currentdir;

system("date > $logfile");
system("hostname >> $logfile");

#-----------------------------------------------#
#  get list of input (nexus) files              #
#-----------------------------------------------#

#my @gene = glob("$seqdir/*");
chdir($seqdir) or die "can't go to sequence directory";
my $genefiles = `ls`;
my @genes = split(/\n/, $genefiles);
my @generoots;
foreach my $gene (@genes){
    my $generoot = $gene;
    $generoot =~ s/\.\w{3}//;
    push @generoots, $generoot;
}
my $nloci = scalar(@genes);
chdir($currentdir) or die "can't go back to original directory";

#-----------------------------------------------#
#  convert nexus to phylip files                #
#  interleaved format: do *not* repeat taxon names
#-----------------------------------------------#

if ($convertphylip) {
  for my $ig (0 .. $#genes){
    my $infn = "$seqdir/${genes[$ig]}";
    my $oufn = "$phylipdir/${generoots[$ig]}.phy";
    my $read = 0;
    my $removeNames = 0; my $nReadNames = 0;
    my $ntax = 0;
    my $nchar = 0;
    open my $FHi, $infn or die "can't open NEXUS gene sequence file";
    open my $FHo, ">", $oufn or die "can't open PHYLIP gene sequence file";
    while (<$FHi>){
	  if ($read){
	    if (/^\s*;/){ last; } # end of alignment
	    if (/^\s*$/){ next; } # ommit blank lines: RAxML doesn't like them
	    if ($removeNames){
	        if (/^[^\s]+\s+(.*)/) { print $FHo "$1\n"; }
	    } else {                    print $FHo $_;
	        $nReadNames++;
	        if ($nReadNames == $ntax){ $removeNames = 1; }
	    }
	  }
	  if ($read==0){
	    if (/ntax\s*=\s*(\d+)/i){
	    	$ntax = $1;
	    	print $FHo " $ntax ";
	    }
	    if (/nchar\s*=\s*(\d+)/i){
		    $nchar = $1;
		    if ($ntax==0){ print "problem in file $infn: found nchar before ntax\n"}
		    print $FHo "$nchar\n";
	    }
	  }
	  if (/^\s*matrix/i){
	    $read=1;
	    if ($ntax==0 or $nchar==0){
		    print "problem in file $infn: was unable to find ntax ($ntax) or nchar ($nchar)\n";
	    }
	  }
    }
    close $FHi;
    close $FHo;
  }
  $seqdir = $phylipdir;
  for my $ig (0..$#genes){
      $genes[$ig] = $generoots[$ig]. ".phy";
  }
}

#-----------------------------------------------#
#  run RAxML on phylip files                    #
#-----------------------------------------------#

open FHlog, ">> $logfile";
chdir($raxmldir) or die ("can't go to raxml directory $raxmldir\n");

for my $ig (0 .. $#genes){
    my $infn = "../$seqdir/${genes[$ig]}";
    my $oufn = "${generoots[$ig]}";
    my $str = "$raxml  -T $numCores -m GTRGAMMA --HKY85  -f a -N $numboot";
    # -f a: rapid bootstrap + search for best ML tree
    # seeds: RAxML will fail with almost no info if a seed is 0
    $str .= " -p " . (int(rand(10000))+1); # seed for randomized stepwise addition
    $str .= " -x " . (int(rand(10000))+1); # seed for rapid bootstrapping
    $str .= " -s $infn -n $oufn";      # input/output file names
    print FHlog "starting RAxML for gene $ig...\n";
    print FHlog "$str\n";
    system($str);
}
chdir($currentdir) or die "can't go back to original directory";
system("date >> $logfile");
close FHlog;

#-----------------------------------------------#
#  restructure output files                     #
#-----------------------------------------------#

# archive phylip files
#`tar -czf ${phylipdir}.tgz $phylipdir`

# delete info files created by raxml
`rm $raxmldir/RAxML_info\.*`;

# delete ".reduced" files created by RAxML for genes with duplicated sequences
my $nreduced=`ls $seqdir | grep ".reduced" | wc -l`; # number of ".reduced" files
$nreduced =~ s/^\s*|\s*$//g; # rm leading & trailing spaces: get the number only
if ($nreduced){
    print "RAxML created $nreduced reduced alignment files: removing them\n";
    `rm $seqdir/*.reduced`;
}

# create directory to contain the bootstrap trees for all genes: 1 file per gene
my $bootpath = "$raxmldir/bootstrap";
make_path $bootpath unless(-d $bootpath);
`mv $raxmldir/RAxML_bootstrap.* $bootpath`;
# do not tar: needed for ASTRAL

# create directory to contain the consensus trees: 2 files per gene
`tar -czf $raxmldir/contrees.tgz $raxmldir/RAxML_bipartitions*`;
`rm -f $raxmldir/RAxML_bipartitions*`;

# create file listing all best trees: one line per gene
my $raxmlOUT = "$raxmldir/besttrees.tre";
`cat $raxmldir/RAxML_bestTree\.* > $raxmlOUT`;
`tar -czf $raxmldir/besttrees.tgz $raxmldir/RAxML_bestTree\.*`;
`rm -f $raxmldir/RAxML_bestTree\.*`;

# ----------------------------------------------#
#   astral analysis                             #
# ----------------------------------------------#

$astraldir = "astral" if !defined($astraldir);
my $bsfile =  "$astraldir/BSlistfiles";
my $astralLOG =  "$astraldir/astral.screenlog";
my $astralOUT =  "$astraldir/astral.tre";

`ls -d $bootpath/* > $bsfile`;

my $astralcmd = "java -jar $astral -i $raxmlOUT -b $bsfile -r $numboot -o $astralOUT > $astralLOG 2>&1";
open FHlog, ">> $logfile";
if ($doastral){
    print FHlog "running astral:\n";
} else {
    print FHlog "astral could be run with:\n";
}
print FHlog "$astralcmd\n";
close FHlog;
if ($doastral){
    system($astralcmd);
}
