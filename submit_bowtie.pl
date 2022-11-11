#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <input fastq dir> <indices prefix> <output dir> <qsub dir> <stat dir> <dataset> <total jobs fn> <required memory (Mb)> <working dir> <fastq version> <number of parallel jobs\n";
	exit 1;
}

my $idir = $ARGV[0];
my $idx_pref = $ARGV[1];
my $odir = $ARGV[2];
my $qsub_dir = $ARGV[3];
my $statdir = $ARGV[4];
my $dataset = $ARGV[5];
my $total_jobfn = $ARGV[6];
my $req_mem = $ARGV[7];
my $wd = $ARGV[8];
my $fastq_ver = $ARGV[9];
my $jobsize = $ARGV[10];

# S      Sanger          --phred33-quals    key/[0-9]
# X      Solexa          --solexa-quals     key/[0-9]
# I      Illumina 1.3+   --solexa1.3-quals  key/[0-9]
# J      Illumina 1.5+   --solexa1.3-quals  key/[0-9]
# L      Illumina 1.8+   --phred33-quals    key [0-9]...

my $fastq_qual_param;
if ($fastq_ver eq "S" or $fastq_ver eq "L") {
	$fastq_qual_param = "--phred33-quals";
} elsif ($fastq_ver eq "X") {
	$fastq_qual_param = "--solexa-quals";
} elsif ($fastq_ver eq "I" or $fastq_ver eq "J") {
	$fastq_qual_param = "--solexa1.3-quals";
} else {
	die "Invalid fastq version ($fastq_ver). Expecting one of: [SXIJL]\n";
}

# separate qsub dir for task
$qsub_dir = $qsub_dir."_bowtie";
system("rm -rf $qsub_dir");
system("mkdir -p $qsub_dir");

system("rm -rf $odir");
system("mkdir -p $odir");

print "$idir\n";
my @files = <$idir/*>;

my $command_fn = $qsub_dir."/commands";
open(OUT,">", $command_fn);
(@files > 0) or die "directory ".$idir." is empty";
foreach my $file (@files) 
{
	my $base_name = fileparse($file);
	my $ofile = $odir."/".$base_name;

	my $batch_command = "bowtie -t -B 1 -a -m 1 --best --strata --chunkmbs 200 $fastq_qual_param $idx_pref $file $ofile\n";
	print OUT $batch_command;
	
}
close OUT;
my $jobname = "bowtie_".$dataset;
my $command = $wd."/R/submit_batch.r $command_fn $odir $qsub_dir $jobname $jobsize $total_jobfn $req_mem $wd";
print "command: $command\n";
system($command) == 0 or die;

$command = $wd."/lscripts/stats_mapping.pl $qsub_dir $statdir";
system($command) == 0 or die;



