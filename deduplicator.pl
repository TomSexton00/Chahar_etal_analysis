#!/usr/bin/env perl

use strict;
use warnings FATAL => qw ( all );
use File::Path;
use File::Basename;
use POSIX;
$|=1;
die "$0 usage: <PAIR DIR> <STATS DIR> <DEDUP DIR> <LAST SOMATIC CHR>" unless ($#ARGV == 3);
my $idir = $ARGV[0];
my $statdir = $ARGV[1];
my $odir = $ARGV[2];
my $last = $ARGV[3];

my $count = 0;
my $unique = 0;
my @chroms = (1..$last,"X");
my $limit = 300000;
my $filecount = 1;

# clean up
system("rm -rf $odir/*") == 0 or die;
system("mkdir -p $odir") == 0 or die;

open OUT,">",$odir."/".$filecount.".dedup" or die $!;
print "Scanning cis reads...\n\n";

my @files = <$idir/*>;
foreach my $chrom(@chroms) {
	print "Chromosome $chrom : ";
	my %pairs;
	my $chromcount;
	my $chromunique;
	foreach my $tmp(@files) {
		open IN,$tmp or die $!;
		my $head = <IN>;
		while (my $tag = <IN>) {
			chomp $tag;
			my @f = split(/\t/,$tag);
			next unless ($f[0] eq $chrom and $f[3] eq $chrom);
			++$count;
			++$chromcount;
			my $low;
			my $high;
			if ($f[1] < $f[4]) {
				$low = $f[1].$f[2];
				$high = $f[4].$f[5];
			}
			else {
				$high = $f[1].$f[2];
				$low = $f[4].$f[5];
			}
			my $id = $low."_".$high;
			next if (exists $pairs{$id});
			++$unique;
			++$chromunique;
			print OUT $tag,"\n";
			$pairs{$id} = 1;
			if ($unique % $limit == 0) {
				close OUT or die $!;
				++$filecount;
				open OUT,">",$odir."/".$filecount.".dedup" or die $!;
			}
		}
	close IN or die $!;
	print ".";
	}
	print 100*($chromunique/$chromcount)," % unique\n";
}

print "\nScanning trans reads : ";
my %pairs;
my $transcount;
my $transunique;
foreach my $tmp(@files) {
	open IN,$tmp or die $!;
	my $head = <IN>;
	while (my $tag = <IN>) {
		chomp $tag;
		my @f = split(/\t/,$tag);
		next unless ($f[0] ne $f[3]);
		++$count;
		++$transcount;
		my $id;
		my @chromcombo = ($f[0],$f[3]);
		@chromcombo = sort {$a cmp $b} @chromcombo;
		if ($chromcombo[0] eq $f[0]) {
			$id = $f[0]."_".$f[1].$f[2].":".$f[3]."_".$f[4].$f[5];
		}
		else {
			$id = $f[3]."_".$f[4].$f[5].":".$f[0]."_".$f[1].$f[2];
		}
		next if (exists $pairs{$id});
		++$unique;
		++$transunique;
		print OUT $tag,"\n";
		if ($unique % $limit == 0) {
			close OUT or die $!;
			++$filecount;
			open OUT,">",$odir."/".$filecount.".dedup" or die $!;
		}
		$pairs{$id} = 1;
	}
	close IN or die $!;
	print ".";
}
print 100*($transunique/$transcount)," % unique\n\n";
close OUT or die $!;

open OUT,">",$statdir."/deduplication.stats" or die $!;
my $per = 100 * ($unique/$count);
print OUT "$unique out of $count are unique - ($per %)\n";
close OUT or die $!;
print "Done!\n";

