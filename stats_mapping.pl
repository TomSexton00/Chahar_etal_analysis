#!/usr/bin/env perl

use strict;
use warnings FATAL => qw ( all );
use File::Path;
use File::Basename;
use POSIX;

die "$0 usage: < MAP_Bowtie_QSUB DIR> <STATS_DIR> " unless ($#ARGV == 1);

my $idir = $ARGV[0];
my $odir = $ARGV[1];
my $total;
my $aligned;
# Store only Stdout-files in the @files array using glob

my @files = glob "$idir/*.stdout";
for (0..$#files){
  $files[$_] =~ s/\.txt$//;
  open (IN,$files[$_]) || die;
  while (<IN>) {
	chomp;
	my $pos= index($_,"# reads processed:"); 
	if ( $pos != -1) {
		my @line= split(/ /);
		$total+=$line[3];
			
	}
	
	my $fuse =index ($_,"Reported");
	if ($fuse != -1) {
		my @line= split(/ /);
		$aligned+=$line[1];
	}


  }
  close (IN);
}
my $ofile=$odir."/mapping.stats";
open (OUT,">",$ofile)|| die;
print OUT "Total reads: ",$total,"\n";
print OUT "Aligned reads: ",$aligned,"\n";
print OUT "Percentage: ",$aligned *100 / $total,"\n";


