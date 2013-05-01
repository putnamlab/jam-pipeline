#!/usr/bin/perl -w

# mm_links2amos.pl
# Convert our read-to-contig hits list based on major-minor+sliced kmers
# to format Amos/Bambus can understand. (Their contigs & reads format,
# but just the headers).
# Our format (note that SNPmers are ignored in producing Bambus format)
# ==========
# lib.ends.IDinLib(<rangeinfo>)	contignum[<crangeinfo>]diag:<S><E>[<rcrangeinfo>]#num	{repeats...}
# OR
# lib.IDinLib.ends[<rangeinfo>](<snpmerpairs>)	contignum[<crangeinfo>]diag:<S><E>[<rcrangeinfo>]#num(<snpmers>)	{repeats...}
# Where the parenthesized <snpmerpairs> and <snpmers> are optional.
#	EXAMPLE, form 1
#		p6.m.000000107(4,46;13,72)      6346883[1,143]47:-f[4,46]#4     1188166[1,69]-12:+r[13,72]#4
#		where
#		lib=	p6	(a sequencing library)
#		ends=	m	(indicating forward and reverse, alternatives are 'f' or 'r' for only one end)
#		IDinLib=	000000107
#		<rangeinfo> => fstart,fend;rstart,rend => 4,46;13,72
#                  (starting locations of first & last interesting kmers in each read;
#                   pay more attention to hit-kmer ranges)
#	                 (if no reverse reads, 
#		contignum=	6346883
#		Crange=	1,143 (starting locations of first & last kmers in contig, whole contig)
#		diag=47	(sum or difference of contig/read offsets)
#		S= -	('-' indicates antidiagonal, so diag field is sum of contig & read locs)
#       (alternatively, '+' would indicate regular diagonal, given as contig_loc - read_loc)
#		E= f (which End of read = forward)
#		Rrange= 4,46 (starting locations of first & last kmers for this read-contig hit, hit portion of read only)
#				(Note that first read hit-kmer corresponds to last contig hit-kmer for antidiagonal hit,
#					so those correspond to kmers starting at 43 (=47-4) and 1 (=47-46), respectively.
#					If we're talking about bases, base 4 of the read is contig base 65 (=47-4+(k-1)) and
#					base 46 of the read is contig base 23 (=47-46+(k-1)).)
# EXAMPLE, form 2
#   10.000000034.m[3,43;3,3](63c4902fac0:63c4982fac0,208defc507f9:208defc507d9)	253261[1,199]47:+f[3,43]#4(63c4902fac0,208defc507f9)	205042[1,438]435:+r[3,3]#1
#   where
#   <snpmerpairs> => 63c4902fac0:63c4982fac0,208defc507f9:208defc507d9
#                    present, as first member of pair, in the read
#   <snpmers> => 63c4902fac0,208defc507f9
#                present as either member of pair, in the contig
#   Note that the omitted parentheses and SNPmer list for the second contig implies an empty list
#   SNPmers are ignored in producing the Bambus output.

#
#	Amos/Bambus format, headers only
# ================================
# Contigs:
# ##contig_name #reads #bases bases, CHKSUM checksum.
# e.g.
# ##asm_10 53 5213 bases, 59DBC21B checksum.
# 
# Reads:
# #read_name(zero_based_offset) [sense] #bases bases, CHKSUM checksum. {rstart rend} <cstart cend>
# e.g.
# #TBCNC59.y1(105) [] 273 bases, ED2069E9 checksum. {1 273} <106 378>
# #TBCNG63.y1(139) [RC] 290 bases, EDFD0386 checksum. {290 1} <140 429>
# Note:
#   * empty sense -> read aligns in forward direction (to top strand of contig)
#     'RC' sense -> read aligns in backward direction (to bottom strand)
#                   and read base numbers are in {rend rstart} order
#   * contig starting base cstart is one greater than zero_based_offset value.
# 
# Conversion:
# ===========
# 1) Generate contig-read lines from our format, with Amos/Bambus format line prefixed by contig info
# contignum #bases read_id(zero_based_offset) [sense] #bases bases, CHKSUM checksum. {rstart rend} <cstart cend>
#    * #reads will be counted from the number of these records
#    * CHKSUM will be zero
#	   * contignum = contignum from our format
#    * #bases = Crange.end + (k-1)
#    * read_id = lib_IDinLib.end (where end = 'f' or 'r') Reformatted from our format because insert name (lib_IDinLib) 
#                needs to be a contiguous string
#    * zero_based_offset = if S (diagonal sense) == '+', (diag + Rrange.start - 1)
#                          else /S == '-'/, (diag - (Rrange.end + 1))
#    * [sense] = if S == '+', []
#                else S == '-', [RC]
#    * #bases = largest known value for this read's Rrange.end + (k-1)
#    * CHKSUM = 0
#    * rstart = if S == '+', Rrange.start
#               else S == '-', Rrange.end + (k-1)
#    * rend = if S == '+', Rrange.end + (k-1)
#             else S == '-', Rrange.start
#    * cstart = zero_based_offset + 1
#    * cend = if S == '+', (diag + Rrange.end) + (k-1)
#             else S == '-', (diag - Rrange.start) + (k-1)

use strict;
use Carp;

use Getopt::Long;
use List::Util;

my $help = 0;
my $k    = 23;
die "Options problem\n" unless
  &GetOptions(
							"k=i"     => \$k,
							"help|h"	=> \$help,
						 );

# Report statistics
#
# By library
# Consistent same-read links 

# sub byDescKcount() {
#  my ($acount) = $a =~ /#(\d+)$/;
#  my ($bcount) = $b =~ /#(\d+)$/;
#  $bcount <=> $acount;
# }

sub check() {
  my $aref = shift;
  # my @list = sort byDescKcount @$aref;
  my @discard = (0) x scalar(@{$aref});
 OUTER: for my $i (0 .. $#{$aref}) {
    my ($contig1, $min1, $max1, $diag1, $strand1, $readend1, $start1, $stop1, $count1) =
      ($$aref[$i] =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)([(]|$)/);
    for my $j (($i+1) .. $#{$aref}) {
      # if $hits overlap, reject lower-scoring or both if too close
      my ($contig2, $min2, $max2, $diag2, $strand2, $readend2, $start2, $stop2, $count2) =
				($$aref[$j] =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)([(]|$)/);
      # Mappings shouldn't overlap in read kmer space
			unless ($stop1 < $start2 or $stop2 < $start1) {
				my ($span1, $span2) = ($stop1 - $start1, $stop2 - $start2); # off-by-1 not important in comparison
				if ($count1 >= $count2 or ($span1 >= $span2)) {
					$discard[$j] = 1;
				}
				elsif ($count2 >= $count1 or ($span2 >= $span1)) {
					$discard[$i] = 1;
				}
			}
      # Contigs tied to read by mappings shouldn't overlap in kmers
      # (kmer-contigs were constructed so as to never share marker kmers)
    }
  }
	my @keep = ();
	for my $ii (0 .. $#{$aref}) {
		push @keep, $$aref[$ii]
			unless $discard[$ii];
	}
  @$aref = @keep;
}

sub amosformat {
	my ($template, $end, $size, $hit) = @_;

	# print "HIT: /$hit/\n";

	my ($contig, $min, $max, $diag, $strand, $readend, $start, $stop, $count) =
		($hit =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)([(]|$)/);
	# print STDERR "contig=$contig, min=$min, max=$max, diag=$diag, strand=$strand, readend=$readend, start=$start, stop=$stop, count=$count\n";
	# print "FIELDS: /$contig/$min/$max/$diag/$strand/$readend/$start/$stop/$count/\n";

	# contignum #contigbases read_id(zero_based_offset) [sense] #readbases bases, CHKSUM checksum. {rstart rend} <cstart cend>
	my $zero_based_offset = ($strand eq '+')? ($diag + $start - 1) : ($diag - ($stop + 1));

	return sprintf("%d %d %s.%s %d [%s] %d bases, 0 checkum. {%d %d} <%d %d>\n",
								 $contig,            # contignum
								 ($max + $k - 1),    # #bases in contig
								 $template, $end,    # for readid = templateid_(f|r)
								 $zero_based_offset,
								 ($strand eq '+') ? "" : "RC", # sense (regular or ReverseComp)
								 $size,              # #bases in read
								 ($strand eq '+')? $start : ($stop + $k - 1), # start of hit in read
								 ($strand eq '+')? ($stop + $k - 1) : $start, # end of hit in read
								 $zero_based_offset + 1,                      # start of hit in contig
								 ($strand eq '+')? ($diag + $stop + $k - 1) : ($diag + ($k - 1) - $start) # end of hit in contig
								);
}

while (<>) {
  chomp;
  my @v = split;
  my $readinfo = shift @v;
	# Can't have a link unless the read or insert has at least two hits
	next unless (@v > 1);

	# Note: ignoring any parenthesized SNPmers pairs at end of readinfo.
	# Note also: sometimes read ranges are set off with parens, sometimes with square brackets
	# print "READINFO: /$readinfo/\n";;
  my ($read, $lib, $ends, $IDinLib, $delim, $range1s, $range1e, $range2s, $range2e) =
		($readinfo =~ /^((\w+)[._](\w+)\.(\w+))(\(|\[)(\d+),(\d+)[];:)](\d*),?(\d*)/);
	# print "FIELDS: /$read/$lib/$ends/$IDinLib/$delim/$range1s/$range1e/$range2s/$range2e/\n";;
	if (length($IDinLib) == 1) {
		# swapped in some variants of this format
		($IDinLib, $ends) = ($ends, $IDinLib);
	}
	# print STDERR "parsed readinfo: read=$read, lib=$lib, ends=$ends, IDinLib=$IDinLib, ranges: $range1s,$range1e,$range2s,$range2e\n";

	my $template  = $lib . "_$IDinLib";
	my ($readFsize, $readRsize);
	if ($ends eq 'm') {
		$readFsize = $range1e + $k - 1;
		$readRsize = $range2e? ($range2e + $k - 1) : 0;
	}
	elsif ($ends eq 'f') {
		$readFsize = $range1e + $k - 1;
		$readRsize = 0;
	}
	else { # 'r'
		$readFsize = 0;
		$readRsize = $range1e + $k - 1;
	}

	my $bestFwdHit    = "";
	my $bestFwdCount  = 0;
	my $bestFwdSpan   = 0;
	my $bestFwdContig = "";
	my $bestRevHit    = "";
	my $bestRevCount  = 0;
	my $bestRevSpan   = 0;
	my $bestRevContig = "";

  for my $hit (@v) {
    my ($contig, $cmin, $cmax, $diag, $strand, $readend, $start, $stop, $count) =
      ($hit =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)([(]|$)/);
    unless (defined $contig and defined $count) {
      print STDERR "PARSEFAIL: $contig $count $hit\n";
      next;
    }
    if ($readend eq 'f') {
      if ($count > $bestFwdCount) {
				$bestFwdCount  = $count;
				$bestFwdSpan   = $stop - $start;
				$bestFwdHit    = &amosformat($template, 'f', $readFsize, $hit);
				$bestFwdContig = $contig;
      }
    }
    else { # 'r'
      if ($count > $bestRevCount) {
				$bestRevCount  = $count;
				$bestRevSpan   = $stop - $start;
				$bestRevHit    = &amosformat($template, 'r', $readRsize, $hit);
				$bestRevContig = $contig;
      }
    }
  }
	print $bestFwdHit, $bestRevHit 
		if ($bestFwdSpan >= $k and $bestRevSpan >= $k and $bestFwdContig ne $bestRevContig);
}

carp "usage: $0 <args>\n" if $help;
