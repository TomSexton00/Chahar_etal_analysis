#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <split_dir> <trunc_dir> <qsub_dir> <junction> <dataset> <total jobs fn> <working dir> <number of parallel jobs> <stats_dir>\n";
	exit 1;
}

my $idir = $ARGV[0];
my $odir = $ARGV[1];
my $qsub_dir = $ARGV[2];
my $junction = $ARGV[3];
my $dataset = $ARGV[4];
my $total_jobfn = $ARGV[5];
my $wd = $ARGV[6];
my $jobsize = $ARGV[7];
my $statdir = $ARGV[8];

# separate qsub dir for task
$qsub_dir = $qsub_dir."_trunc";
system("rm -rf $qsub_dir");
system("mkdir -p $qsub_dir");

system("rm -rf $odir");
system("mkdir -p $odir");

system("rm -rf $statdir");
system("mkdir -p $statdir");

print "$idir\n";
my @files = <$idir/*>;

my $command_fn = $qsub_dir."/commands";
open(OUT,">", $command_fn);
(@files > 0) or die "directory ".$idir." is empty";
foreach my $file (@files) 
{
	my $base_name = fileparse($file);
	my $ofile = $odir."/".$base_name;

	my $batch_command = "$wd/lscripts/trunc.pl $file $ofile $junction";
	print OUT $batch_command,"\n";
	
}
close OUT;
my $jobname = "trunc_".$dataset;
my $command = $wd."/R/submit_batch.r $command_fn $odir $qsub_dir $jobname $jobsize $total_jobfn 0 $wd";
print "command: $command\n";
system($command) == 0 or die;

$command = $wd."/lscripts/stats_truncation.pl $qsub_dir $statdir";
system($command) == 0 or die $!;
