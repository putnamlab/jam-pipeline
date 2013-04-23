#!/usr/bin/perl -w

# fastq2fam.pl
# Reformats FASTQ format (insisting on exactly 4 lines per read) to 
# masked FASTA: low-quality bases converted
#   * to 'N' (if < $minqual)
#   * or lowercase (if >= $minqual but < $softmask)
#
# Does not filter or reject reads, just deprecates or hides the 
# low-quality parts.
#
# Note also the $zeroqual option, to deal with different versions of
# FASTQ format using a different ASCII character for the zero-quality 
# start of the quality scale.
#   * typically 33 (ASCII '!') or 64 (ASCII '@')
#   * see http://en.wikipedia.org/wiki/FASTQ_format
#
# Other options relate to renaming reads with:
#   * a required pattern $readpat for sanity-checking of read headers 
#     (you will want to customize)
#   * a uniform $prefix for this invocation of script 
#     (e.g. a library name and/or run identifier)
#   * batching parameters $batch and $bsize, used for numbering 
#     reads from multiple files in one contiguous sequence; this batch
#     would (if full) get read numbers
#         ($batch - 1)*$bsize .. ($batch*$bsize - 1)
#   * a uniform $suffix for this invocation of script
#     (e.g. ".f" for forward reads, ".r" for reverses)
#
# Reads from standard input and writes to standard output for maximum flexibility
# of what the original and new files are named.

use strict;
use Carp;

use Getopt::Long;

my $help     = 0;
my $minqual  = 20;
my $softmask = 30;
my $zeroqual = 33; # as for original FASTQ and Illumina pipeline 1.8+
my $prefix   = "";
my $suffix   = "";
my $batch    = 1;
my $bsize    = 4e6;
my $readpat  = "HWI|FC[CD]";

die "Options problem\n" unless
  &GetOptions(
              "help|h"          => \$help,
              "minqual|m=i"     => \$minqual,
              "softmask|soft=i" => \$softmask,
              "zeroqual|zq=i"   => \$zeroqual, # ASCII code of zero quality value
              "prefix|pre|p=s"  => \$prefix,  # prefix for all reads, e.g. common library
              "suffix|suf=s"    => \$suffix,  # suffix for all reads, e.g. .f or .r
              "readpattern|r=s" => \$readpat, # required read header pattern
              "batch|bnum=i"    => \$batch,   # batch number for read file within library (default 4e6 reads each)
              "batchsize|bs=i"  => \$bsize,   # batch size (multiply by $batch - 1 to get $readbase)
             );

my $minq  = chr($zeroqual + $minqual);
my $softm = chr($zeroqual + $softmask);
carp "usage: $0 <args>\n" if $help;

# Line and read counters
my ($line, $nreads) = (0, 0);
# Hold the lines for descriptions or DNA bases
my ($descrip, $bases) = ("", "");

my $readbase = ($batch - 1)*$bsize;
while (<>) {
  $line++;
  my $part = $line % 4;
  if (1 == $part) {
    # First line of four is description/header of the read.
    die unless /^@($readpat)/;
    $nreads++;
    if ($prefix or $readbase or $suffix) {
      my $readid = sprintf("%s.%09d.%s", $prefix, $readbase + $nreads - 1, $suffix);
      s/^@/>$readid /;
    }
    else {
      s/^@/>/; 
    }
    print; # the description line as modified by the s/// expressions above.
    $descrip = $_;
  } 
  else { 
    if (2 == $part) {
      # Second line of four is the DNA bases
      $bases = $_;
    }
    elsif (0 == $part) {
      # Fourth line of four is ASCII-encoded quality scores
      chomp;
      my $i = 0;
      for my $q (split //) {
        # print "$i: q$q for", substr($bases, $i, 1), "\n";
        if ($q lt $minq) {
          substr($bases, $i, 1, 'N');
        }
        elsif ($q lt $softm) {
          substr($bases, $i, 1, lc(substr($bases, $i, 1)));
        }
        else {
          # In case input sequences were all lower-case
          substr($bases, $i, 1, uc(substr($bases, $i, 1)));
        }
        $i++;
      }
      # The quality-masked DNA string
      print $bases;
    }
    # third line of four is a description for the quality scores, usually left blank

    # Uncomment this to get a progress report.
    # print STDERR "Nreads: $nreads\r" unless $nreads%1e5;
  } 
}
