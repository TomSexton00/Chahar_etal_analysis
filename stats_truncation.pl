#!/usr/bin/env perl

use strict;
use warnings FATAL => qw ( all );
use File::Path;
use File::Basename;
use POSIX;

die "$0 usage: < MAP_Trunc_QSUB DIR> <stats dir> " unless ($#ARGV == 1);

my $idir = $ARGV[0];
my $odir = $ARGV[1];
my $total;
my $truncated;
my $truncation_length;
my $compteur_files;
# Store only Stdout-files in the @files array using glob

my @files = glob "$idir/*.stdout";
for (0..$#files){
  $files[$_] =~ s/\.txt$//;
  open (IN,$files[$_]) || die;
  while (<IN>) {
	chomp;
	my $pos= index($_,"reads truncated"); 
	if ( $pos != -1) {
		my @line= split(/ /);
		$total+=$line[3];
		$truncated+=$line[0];
			
	}
	
	my $fuse =index ($_,"Average truncation length:");
	if ($fuse != -1) {
		my @line= split(/ /);
		$truncation_length+=$line[3];
	}


  }
  $compteur_files++;
  close (IN);
}
my $ofile=$odir."/truncation.stats";
open (OUT,">",$ofile)|| die;
print OUT "Total reads: ",$total,"\n";
print OUT "Truncated reads: ",$truncated,"\n";
print OUT "Percentage: ",$truncated *100 / $total,"\n";
print OUT "Average truncation length: ",$truncation_length/$compteur_files,"bp","\n";

close (OUT);
