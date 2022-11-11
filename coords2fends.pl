#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use Data::Dumper;

if ($#ARGV == -1) {
	print STDERR "usage: $0 <fends> <read length> <discard facing threshold> <max segment length> <output mat prefix> <input file1> [input file2] ...\n";
	exit 1;
}

my $in_fends_fn = $ARGV[0];
my $read_len = $ARGV[1];
my $discard_facing_threshold = $ARGV[2];
my $max_len = $ARGV[3];
my $mat_fn_prefix = $ARGV[4];
shift;shift;shift;shift;shift;
my @ifns =  @ARGV;
print STDERR "Input files: ", join("\n", @ifns), "\n";

# table with fends
our %fends;

# hash table to quickly get from approximate coord to fend
our %coord2index;

##########################################################################################
# read fends file
##########################################################################################

open(IN, $in_fends_fn) || die $in_fends_fn;
print STDERR "Reading input file $in_fends_fn into hash...\n";
my $header = <IN>;
my %h = parse_header($header);
while (my $line = <IN>) {
	chomp $line;
	my @f = split("\t", $line);
	my $fend = $f[$h{fend}];
	my $chr = $f[$h{chr}];
	my $coord = $f[$h{coord}];
	my $strand = $f[$h{strand}];
	my $frag = $f[$h{frag}];
	my $frag_len = $f[$h{frag_len}];

	!defined($fends{$fend}) or die "non-unique fend";
	$fends{$fend} = {};
	$fends{$fend}->{fend} = $fend;
	$fends{$fend}->{frag} = $frag;
	$fends{$fend}->{strand} = $strand;
	$fends{$fend}->{chr} = $chr;
	$fends{$fend}->{coord} = $coord;
	$fends{$fend}->{frag_len} = $frag_len;

	if (!defined($coord2index{$chr}))
	{
		$coord2index{$chr} = {};
		$coord2index{$chr}->{coords} = {};
	}
	$coord2index{$chr}->{coords}->{$coord} = $fend;
}
close(IN);

# compute sorted coords per chrom
for my $chr (keys %coord2index)
{
	my @sorted = sort {$a <=> $b} keys %{$coord2index{$chr}->{coords}};
	$coord2index{$chr}->{sorted_coords} = \@sorted;
}

##########################################################################################
# parse pair files
##########################################################################################

my $read_count = 0;
my $discard_reads = 0; # reads that cannot be typed (do not have two valid fends, e.g. are located on end of chr)
my $one_side_cutter = 0;

my @types = ("s0", "s1", "s2");
my %fend_matrix;
my %read_stats;

foreach my $type (@types) {
	# out containers
	$fend_matrix{$type} = {};
	# read stats
	$read_stats{$type} = {};
	$read_stats{$type}->{count} = 0;
	$read_stats{$type}->{no_ligation} = 0;
	$read_stats{$type}->{self_ligation} = 0;
	$read_stats{$type}->{no_restriction} = 0;
	$read_stats{$type}->{seglen_out_of_range} = 0;
	$read_stats{$type}->{cis_0_1k} = 0;
	$read_stats{$type}->{cis_1k_10k} = 0;
	$read_stats{$type}->{cis_10k_100k} = 0;
	$read_stats{$type}->{cis_100k_1m} = 0;
	$read_stats{$type}->{cis_1m_max} = 0;
	$read_stats{$type}->{trans} = 0;
}

foreach my $ifn (@ifns) {

	print STDERR "traversing file: $ifn\n";
	open(IN, $ifn) || die $ifn;

	$header = <IN>;
	%h = parse_header($header);

	# if header is not defined we add it here
	if (!defined($h{coord1})) {
		$h{chr1} = 0;
		$h{coord1} = 1;
		$h{strand1} = 2;
		$h{chr2} = 3;
		$h{coord2} = 4;
		$h{strand2} = 5;
	}

	while (my $line = <IN>) {
		$read_count++;
		print STDERR "line: $read_count\n" if ($read_count % 100000 == 0);
		chomp $line;
		my @f = split("\t", $line);

		my $chr1 = $f[$h{chr1}];
		my $coord1 = $f[$h{coord1}];
		my $strand1 = $f[$h{strand1}];
		my $chr2 = $f[$h{chr2}];
		my $coord2 = $f[$h{coord2}];
		my $strand2 = $f[$h{strand2}];

		# skip undefined chroms
		next if !defined($coord2index{$chr1}) || !defined($coord2index{$chr2});

		# trim read by 4 bps, i.e., move coord to point to the first nt after the (expected) cutter site
		$coord1 += 4 if ($strand1 eq "+");
		$coord2 += 4 if ($strand2 eq "+");
		$coord1 += $read_len-5 if ($strand1 eq "-");
		$coord2 += $read_len-5 if ($strand2 eq "-");

		my ($fend1, $fwd_dist1, $bck_dist1) = coord_to_fend($chr1, $coord1, $strand1);
		my ($fend2, $fwd_dist2, $bck_dist2) = coord_to_fend($chr2, $coord2, $strand2);
		if ($fend1 == -1 or $fend2 == -1)
		{
			$discard_reads++;
			next;
		}

		# identify type of restriction
		my $type;
		my $on_cutter_site1 = ($bck_dist1 == 2);
		my $on_cutter_site2 = ($bck_dist2 == 2);
		if ($on_cutter_site1 && $on_cutter_site2) {
			$type = "s2";
		} elsif (!$on_cutter_site1 && !$on_cutter_site2) {
			$type = "s0";
		} else {
			$type = "s1";
		}

		# to estimate chance of one side to fall on cutter
		$one_side_cutter++ if ($on_cutter_site1);

		# total read count for type
		$read_stats{$type}->{count}++;

		# compute stats for different events
		my $dist = abs($fends{$fend1}->{coord} - $fends{$fend2}->{coord});
		my $frag1 = $fends{$fend1}->{frag};
		my $frag2 = $fends{$fend2}->{frag};
		my $face_towards = (($coord1 < $coord2) && ($strand1 eq "+") && ($strand2 eq "-")) ||
			(($coord2 < $coord1) && ($strand1 eq "-") && ($strand2 eq "+"));

		# skip if ligation is within a single fragment 
		if ($frag1 == $frag2) {
			if ($face_towards) {
				$read_stats{$type}->{no_ligation}++;
			} else {
				$read_stats{$type}->{self_ligation}++;
			}
			next;
		} 

		# also, discard the close cis facing products - they are probably no restriction events
		if (($chr1 eq $chr2) && $face_towards && ($dist < $discard_facing_threshold)) {
			$read_stats{$type}->{no_restriction}++;
			next;
		}

		# finally if inferred segment distance is too long we discard read as well
        # this is mainly relevant for 6-cutters
		my $segment_len = $fwd_dist1 + $fwd_dist2;
		if ($max_len > 0 && $segment_len > $max_len)
		{
			$read_stats{$type}->{seglen_out_of_range}++;
			next;
		}

		# now, breakdown the normal ligation events according to distance
		if ($chr1 eq $chr2) {
			if ($dist < 1000) {
				$read_stats{$type}->{cis_0_1k}++;
			} elsif ($dist < 10000) {
				$read_stats{$type}->{cis_1k_10k}++;
			} elsif ($dist < 100000) {
				$read_stats{$type}->{cis_10k_100k}++;
			} elsif ($dist < 1000000) {
				$read_stats{$type}->{cis_100k_1m}++;
			} else {
				$read_stats{$type}->{cis_1m_max}++;
			}
		} else {
			$read_stats{$type}->{trans}++;
		}

		# track number of covered fends
		my $fend_small = ($fend1 <= $fend2) ? $fend1 : $fend2;
		my $fend_large = ($fend1 <= $fend2) ? $fend2 : $fend1;
		$fend_matrix{$type}->{$fend_small} = {} if !defined($fend_matrix{$type}->{$fend_small});
		$fend_matrix{$type}->{$fend_small}->{$fend_large} = 0 if !defined($fend_matrix{$type}->{$fend_small}->{$fend_large});
		$fend_matrix{$type}->{$fend_small}->{$fend_large}++;
	}
	close(IN);
}

######################################################################################################
# write mat file
######################################################################################################

foreach my $type (@types) 
{
	my $mat_fn = $mat_fn_prefix."_$type.mat";
	print STDERR "writing file: $mat_fn\n";
	open(OUT, ">", $mat_fn) || die;
	print OUT "fend1\tfend2\tcount\n";
	foreach my $fend1 (sort { $a <=> $b } keys %{$fend_matrix{$type}})
	{
		foreach my $fend2 (sort { $a <=> $b } keys %{$fend_matrix{$type}->{$fend1}})
		{
			my $count = $fend_matrix{$type}->{$fend1}->{$fend2};
			print OUT $fend1, "\t" ,$fend2, "\t", $count, "\n";
		}
	}
	close(OUT);
}

######################################################################################################
# Write stats
######################################################################################################

# global stats for all types
my $gstats_fn = $mat_fn_prefix."_all.mat.stats";
open(OUT, ">", $gstats_fn) || die;
print OUT "total_reads\tside1_on_cutter_site\tdiscarded_reads\n";
print OUT $read_count, "\t", $one_side_cutter, "\t", $discard_reads, "\n";
close(OUT);

# compute expected sizes of datasets
my $p = $one_side_cutter / $read_count;
my %expected;
$expected{"s2"} = $p * $p * $read_count;
$expected{"s0"} = (1-$p) * (1-$p) * $read_count;
$expected{"s1"} = 2 * (1-$p) * $p * $read_count;

foreach my $type (@types) 
{
	my $read_stats_fn = $mat_fn_prefix."_$type.read.stats";
	open(OUT, ">", $read_stats_fn) || die;
	print OUT "total\texpected_for_type\tself_ligation\tno_ligation\tno_restriction\tseglength_out_of_range\tcis_0_1k\tcis_1k_10k\tcis_10k_100k\tcis_100k_1m\tcis_1m_max\ttrans\n";
	print OUT $read_stats{$type}->{count}, "\t", int($expected{$type}), "\t";
	print OUT $read_stats{$type}->{self_ligation}, "\t", $read_stats{$type}->{no_ligation}, "\t", $read_stats{$type}->{no_restriction}, "\t", $read_stats{$type}->{seglen_out_of_range}, "\t";
	print OUT $read_stats{$type}->{cis_0_1k}, "\t", $read_stats{$type}->{cis_1k_10k}, "\t", $read_stats{$type}->{cis_10k_100k}, "\t", $read_stats{$type}->{cis_100k_1m}, "\t", 
	$read_stats{$type}->{cis_1m_max}, "\t", $read_stats{$type}->{trans}, "\n";
	close(OUT);
}

print STDERR "coord2fends done\n";

######################################################################################################
# Subroutines
######################################################################################################

# check if C is between and A,B
sub between
{
	my ($A, $B, $C) = @_;
	$B = $B - $A;
	$C = $C - $A;
	return (($B>$C && $C>0) || ($B<$C && $C<0));
}

sub add_fend
{
	my ($chr, $bin, $fend, $strand) = @_;
	$coord2index{$chr} = {} if !defined($coord2index{$chr});
	$coord2index{$chr}->{$bin} = {} if !defined($coord2index{$chr}->{$bin});

	# mark bin if multiple fends map to it
	if (!defined($coord2index{$chr}->{$bin}->{$strand}))
	{	
		$coord2index{$chr}->{$bin}->{$strand} = $fend;
	}
	else 
	{
		$coord2index{$chr}->{$bin}->{$strand} = -1;
	}
}

sub apprx_lines
{
	my ($fn) = @_;
	my $tmp = "/tmp/".getlogin()."_apprx_lines.tmp";
	system("head -n 100000 $fn > $tmp");
	my $size_head = -s $tmp;
	my $size_all = -s $fn;
	$size_head > 0 or die;
	return (int($size_all/$size_head*100000));
}

sub perc_str
{
	my ($n, $total) = @_;
    return ($n." (".(int(1000 * $n / $total)/10)."%)");
}

sub perc_str2
{
	my ($n, $total) = @_;
    return ((int(1000 * $n / $total)/10)."%");
}

# returns first element above/below value in sorted array
sub binary_search 
{
	my $arr = shift;
	my $value = shift;
	my $above = shift; 

	my $left = 0;
	my $right = $#$arr;
	
	while ($left <= $right) {
		my $mid = ($right + $left) >> 1;
		my $c = $arr->[$mid] <=> $value;
		return $mid if ($c == 0);
		if ($c > 0) {
			$right = $mid - 1;
		} else {
			$left  = $mid + 1;
		}
	}
	$left = -1 if ($left > $#$arr);
	$right = -1 if ($right < 0);
	return (($above eq "+") ? $left : $right);
}

sub parse_header
{
	my ($header) = @_;
	chomp($header);
	my @f = split("\t", $header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}

# returns (fend, fwd distance, bck distance)
sub coord_to_fend
{
	my ($chr, $coord, $strand) = @_;
	defined($coord2index{$chr}) or die;

	# opposite strand
	my $ostrand = ($strand eq "+") ? "-" : "+";

	# search fwd and bck in coord list
	my $index_p = binary_search($coord2index{$chr}->{sorted_coords}, $coord, $strand);
	my $index_m = binary_search($coord2index{$chr}->{sorted_coords}, $coord, $ostrand);
	return ((-1,-1,-1)) if (($index_p == -1) || ($index_m == -1));

	# get fwd and bck coords
	my $fcoord_p = $coord2index{$chr}->{sorted_coords}[$index_p];
	my $fcoord_m = $coord2index{$chr}->{sorted_coords}[$index_m];

	# get fwd and bck distance
	my $dist_p = abs($fcoord_p - $coord);
	my $dist_m = abs($fcoord_m - $coord);

	# get fend on fwd side
	defined ($coord2index{$chr}->{coords}->{$fcoord_p}) and defined($fcoord_p) or die;
	my $fend = $coord2index{$chr}->{coords}->{$fcoord_p};

	return ($fend, $dist_p, $dist_m);
}
