#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);

if ($#ARGV == -1) {
	print "usage: $0 <ifn> <ofn> <junction>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $ofn = $ARGV[1];
my $sequence = $ARGV[2];
my $reads;
my $truncated_v;
my @truncated_lengths;

my $sequence_no_underscore = $sequence;
   $sequence_no_underscore =~ s/_//g;    #Remove the '_' delimit symbol from the sequence

#Tom changelog - this way keeps both HindIII sites in the truncation. Should only keep the first one
#my ($length_same) = split( /_/, $sequence );
    #$length_same = length($length_same);    #Number of bases in which the ligation junction and reference genome correspond before first mismatch
my @add = split(/_/,$sequence);


open IN,$ifn or die $!;
open OUT,">",$ofn or die $!;

while (<IN>) {
	if (/^@/ or /^[ATCG]+_/ or /^NO_BARCODE_/) {
		my $line1 = $_;
		my $line2 = scalar <IN>;
		my $line3 = scalar <IN>;
		my $line4 = scalar <IN>;
		chomp $line2;
		chomp $line4;
		++$reads;
  		my $pos = index($line2,$sequence_no_underscore,0);
		if ($pos==-1) {
			print OUT $line1,$line2,"\n",$line3,$line4,"\n";
		}
		else {
	
               		 #Finds the ligation site, if present, and truncates accordingly
                	#(swapping the ligation sequence for the restriction site up until the cut site).
                	#Counts the sequences truncated/not-truncated.
                	my @truncated = split(/$sequence_no_underscore/,$line2);
                    	$line2 = $truncated[0].$add[0];
                    	$line4 = substr( $line4, 0, ( length $line2 ) );
                    	$truncated_lengths[0] += length $line2;
                    	$truncated_lengths[1]++;
                    	$truncated_v++;
		 	print OUT $line1 . $line2 . "\n" . $line3 . $line4 . "\n";
		}
	}
}

close IN or die $!;
close OUT or die $!;
print "$truncated_v out of $reads reads truncated\n";
my $avg = $truncated_lengths[0]/$truncated_lengths[1];
print "Average truncation length: $avg bp\n";

