#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <input dir> <output dir> <omit edge key sequence T|F> <qsub script dir> <edge pref> <dataset> <total jobs fn> <working dir> <fastq version>\n";
	exit 1;
}
my $idir = $ARGV[0];
my $odir = $ARGV[1];
my $omit = $ARGV[2];
my $qsub_dir = $ARGV[3];
my $edge_pref = $ARGV[4] eq "NA" ? "" : $ARGV[4];
my $dataset = $ARGV[5];
my $total_jobfn = $ARGV[6];
my $wd = $ARGV[7];
my $fastq_version = $ARGV[8];

# separate qsub dir for task
$qsub_dir = $qsub_dir."_unite_pairs";
system("rm -rf $qsub_dir");
system("mkdir -p $qsub_dir");

# clean up
system("rm -rf $odir/*") == 0 or die;
system("mkdir -p $odir") == 0 or die;

open(IN, "$idir/files") or die "missing file $idir/files";
my @files = <IN>;
close IN;

my $command_fn = $qsub_dir."/commands";
open(OUT,">", $command_fn);

foreach my $file_pair (@files)
{
	chomp($file_pair);
	my @f = split(" ", $file_pair);
	my $f1 = "$idir/$f[0]";
	my $f2 = "$idir/$f[1]";
	
	my $ofile = $f[0];
	#$ofile =~ /(.*)${edge_pref}[0-9]\.fastq(.*)/;
	$ofile =~ /(.*)${edge_pref}[0-9]\.(.*)/;
	$ofile = $odir."/".$1.".raw".$2;
	my $command = "$wd/lscripts/bowtie_prepare_pairs.pl $f1 $f2 $ofile T $omit $fastq_version";

	print OUT $command, "\n";

}
close OUT;

my $jobname = "unite_".$dataset;
my $jobsize = 50;
my $command = $wd."/R/submit_batch.r $command_fn $odir $qsub_dir $jobname $jobsize $total_jobfn 0 $wd";
print "command: $command\n";
system($command) == 0 or die;
