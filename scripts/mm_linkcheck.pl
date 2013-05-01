#!/usr/bin/perl -w

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
      ($$aref[$i] =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)(\S*)$/);
    for my $j (($i+1) .. $#{$aref}) {
      # if $hits overlap, reject lower-scoring or both if too close
      my ($contig2, $min2, $max2, $diag2, $strand2, $readend2, $start2, $stop2, $count2) =
				($$aref[$j] =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)(\S*)$/);
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


while (<>) {
  my $bestFwdHit = "";
  my $bestFwdCount = 0;
  my $bestFwdSpan  = 0;
  my $otherFwdSum  = 0;
  my $bestRevHit = "";
  my $bestRevCount = 0;
  my $bestRevSpan  = 0;
  my $otherRevSum  = 0;
  my @fkeep = ();
  my @rkeep = ();

  chomp;
  my @v = split;
  my $readinfo = shift @v;
  my ($read, $lib) = ($readinfo =~ /^((\w+)\.(\w+))/);

  @v = sort(@v);

  my $prev = 0;
  for my $hit (@v) {
    my ($contig, $cmin, $cmax, $diag, $strand, $readend, $start, $stop, $count) =
   	  ($hit =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)(\S*)$/);
    unless (defined $contig and defined $count) {
      print STDERR "PARSEFAIL: $contig $count $hit\n";
      next;
    }
    if ($contig != $prev) {
      push @fkeep, $bestFwdHit 
				if ($bestFwdSpan >= 10 and $bestFwdCount >= $otherFwdSum*2 + 1);
      push @rkeep, $bestRevHit 
				if ($bestRevSpan >= 10 and $bestRevCount >= $otherRevSum*2 + 1);
      $bestFwdHit = "";
      $bestFwdCount = 0;
      $bestFwdSpan  = 0;
      $otherFwdSum  = 0;
      $bestRevHit = "";
      $bestRevCount = 0;
      $bestRevSpan  = 0;
      $otherRevSum  = 0;
      $prev = $contig;
    }
    if ($readend eq 'f') {
      if ($count > $bestFwdCount) {
				$otherFwdSum += $bestFwdCount;
				$bestFwdCount = $count;
				$bestFwdSpan  = $stop - $start;
				$bestFwdHit   = $hit;
      }
      else {
				$otherFwdSum += $count;
			}
    }
    else { # 'r'
      if ($count > $bestRevCount) {
				$otherRevSum += $bestRevCount;
				$bestRevCount = $count;
				$bestRevSpan  = $stop - $start;
				$bestRevHit   = $hit;
      }
      else {
				$otherRevSum += $count;
      }
    }
  }
  # Handle last contig hit(s)...
  push @fkeep, $bestFwdHit 
    if ($bestFwdSpan >= 10 and $bestFwdCount >= $otherFwdSum*2 + 1);
  push @rkeep, $bestRevHit 
    if ($bestRevSpan >= 10 and $bestRevCount >= $otherRevSum*2 + 1);
	
  # This is where to check consistency of best hits per contig 
  # (that they don't imply overlaps of the contigs, unless it's just a forward
  #  and reverse pair having shared kmers implying overlap of the same contig)
  &check(\@fkeep);
  &check(\@rkeep);
  my @keep = (@fkeep, @rkeep);

	sub rc {
		my $string = reverse(shift);
		$string =~ y/acgtACGT+-/tgcaTGCA\-+/;
		return $string;
	}

  my $minspace = 500 - ($k - 1);
  # Process hits
  if (@keep > 1) {
    
    for (my $i = 0; $i < @keep; $i++) {
      my ($contig1, $min1, $max1, $diag1, $strand1, $readend1, $start1, $stop1, $count1) =
			  ($keep[$i] =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)(\S*)$/);
      # Space from read start to contig end is
      #    +strand (regular diagonal = contigpos - readpos):
      #       contig_max + 1 - contig_offset_for_read_base1 
      #       = contig_max + 1 - (diag + read_base1)
      #       = contig_max + 1 - (diag + 1)
      #       = contig_max - diag
      #    -strand (antidiagonal = contigpos + readpos)
      #       contig_offset_for_read_base1 + 1 - contig_min
      #       = (diag - read_base1) + 1 - contig_min
      #       = (diag - 1) + 1 - contig_min
      #       = diag - contig_min
      my $tailspace1 = ($strand1 eq '+') ? ($max1 - $diag1) : ($diag1 - $min1);
      for (my $j = $i+1; $j < @keep; $j++) {
				my ($contig2, $min2, $max2, $diag2, $strand2, $readend2, $start2, $stop2, $count2) =
				  ($keep[$j] =~ /^(\d+)\[(\d+),(\d+)\](-?\d+):([+-])([fr])\[(\d+),(\d+)\]#(\d+)(\S*)$/);
				
				my $tailspace2 = ($strand2 eq '+') ? ($max2 - $diag2) : ($diag2 - $min2);
				
				if ($contig1 eq $contig2) {
					# Must be forward & reverse mapping to same contig. Report implied template size,
					# along with flanking sequence lengths,
					# to STDERR only.
					if ($readend1 eq $readend2) {
						# This case should be eliminated by previous stages...
						print STDERR "FAULT: $readinfo $keep[$i] $keep[$j]\n";
						next;
					}
					elsif ($strand1 eq $strand2) {
						print STDERR "O&O_FAIL: $readinfo $keep[$i] $keep[$j]\n";
						next;
					}
					else {
						my $insertsize = 0;
						if ($strand1 eq '+') {
							# and strand2 eq '-'
							$insertsize = ($k - 1) + ($diag2 - $diag1) - 1;
						}
						else { # strand1 eq '-'
							# and strand2 = '+'
							$insertsize = ($k - 1) + ($diag1 - $diag2) - 1;
						}
						if ($insertsize < 101) {
							print STDERR "TOO_CLOSE:\t$insertsize\t$tailspace1\t$tailspace2\t$readinfo\t$keep[$i]\t$keep[$j]\n";
						}
						elsif ($tailspace1 >= $minspace and $tailspace2 >= $minspace) {
							print STDERR "LOOSE_SPACE:\t$insertsize\t$tailspace1\t$tailspace2\t$readinfo\t$keep[$i]\t$keep[$j]\n";
						}
						else {
							print STDERR "TIGHT_SPACE:\t$insertsize\t$tailspace1\t$tailspace2\t$readinfo\t$keep[$i]\t$keep[$j]\n";
						}
					}
				}
				else { # linking different contigs
					my $type;
					
					if ($readend1 eq $readend2) {
						$type = uc($readend1);
						my ($pt1, $pt2);
						if ($stop1 < $start2) {
							if ($strand1 eq '+') {
								$pt1 = $max1 - $diag1;
							}
							else { # '-'
								$pt1 = $diag1 - $min1;
							}
							if ($strand2 eq '+') {
								$pt2 = $min2 - $diag2;
							}
							else { # '-'
								$pt2 = $diag2 - $max2;
							}
							my $gap = $pt2 - $pt1;
							print STDOUT "$contig1\[$min1,$max1\]\t$contig2\[$min2,$max2\]", 
								"\t$type$strand1$strand2\t", $gap,
									"\t$readinfo",
										"\t$start1,$stop1#$count1",
											"\t$start2,$stop2#$count2\n";
							print STDOUT "$contig2\[$min2,$max2\]\t$contig1\[$min1,$max1\]", 
								"\t$type", &rc("$strand1$strand2"), "\t", $gap,
									"\t$readinfo",
										"\t$start2,$stop2#$count2",
											"\t$start1,$stop1#$count1\n";
						}
						elsif ($stop2 < $start1) {
							if ($strand2 eq '+') {
								$pt2 = $max2 - $diag2;
							}
							else { # '-'
								$pt2 = $diag2 - $min2;
							}
							if ($strand1 eq '+') {
								$pt1 = $min1 - $diag1;
							}
							else { # '-'
								$pt1 = $diag1 - $max1;
							}
							my $gap = $pt1 - $pt2;
							print STDOUT "$contig2\[$min2,$max2\]\t$contig1\[$min1,$max1\]", 
								"\t$type$strand2$strand1\t", $gap,
									"\t$readinfo",
										"\t$start2,$stop2#$count2",
											"\t$start1,$stop1#$count1\n";
							print STDOUT "$contig1\[$min1,$max1\]\t$contig2\[$min2,$max2\]", 
								"\t$type", &rc("$strand2$strand1"), "\t", $gap,
									"\t$readinfo",
										"\t$start1,$stop1#$count1",
											"\t$start2,$stop2#$count2\n";
						}
						else {
							croak "Unorderable hits $keep[$i] and $keep[$j]\n";
						}
					}
					else {
						$type = ($readend1 eq 'f')? 'M' : 'm';
						print STDOUT "$contig1\[$min1,$max1\]\t$contig2\[$min2,$max2\]", 
							"\t$type$strand1$strand2\t", $tailspace1 + $tailspace2 + ($k - 1), 
								"\t$readinfo",
									"\t$start1,$stop1#$count1",
										"\t$start2,$stop2#$count2\n";
						$type = ($readend2 eq 'f')? 'M' : 'm';
						print STDOUT "$contig2\[$min2,$max2\]\t$contig1\[$min1,$max1\]",
							"\t$type$strand2$strand1\t", $tailspace1 + $tailspace2 + ($k - 1),
								"\t$readinfo",
									"\t$start2,$stop2#$count2",
										"\t$start1,$stop1#$count1\n";
					}
				}
      }
    }
  }
}

carp "usage: $0 <args>\n" if $help;
