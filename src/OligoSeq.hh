///////////////////////////////////////////////////////////////////////////////
// OligoTools
// = OligoTools for genome analysis
// by Paul Havlak
// copyright 2009-2013 Rice University
//
// OligoTools for genome analysis by Paul Havlak is licensed under a 
// <a rel="license" href="http://creativecommons.org/licenses/by/3.0/deed.en_US">
// Creative Commons Attribution 3.0 Unported License.
// </a>
//
// OligoTools for genome analysis incorporates unpublished modules for Kmer 
// hash tables originally developed under the Atlas Project 
// (Atlas Whole Genome Assembly Suite) at Baylor College of Medicine
// Human Genome Sequencing Center during 2003-2006.
// The BCM-HGSC copyrighted predecessor is available through the 
// Atlas WGA Suite source distribution on the HGSC web site:
//    http://www.hgsc.bcm.tmc.edu/content/atlas-whole-genome-assembly-suite
// 
// Extensive modifications and additional modules, including all features
// for SNPmer pairing and library bitvectors, were developed in the Putnam
// Lab at Rice University:
//    http://nputnam.web.rice.edu
// 
// Open-source repository: part of the Putnam Lab JAM project
//    https://github.com/putnamlab/jam-pipeline
//    https://github.com/putnamlab/jam-pipeline/tree/master/source
//     
// Contact: Paul Havlak:
//    havlak@rice.edu
//    havlak@alumni.rice.edu
//    http://www.linkedin.com/in/havlak
///////////////////////////////////////////////////////////////////////////////
// OligoSeq.hh
// Paul Havlak, March 2004
// $Header$
// Provide a stream of oligos (k-mers for a constant k)
// read from a DNA sequence file.  Each location in the DNA
// sequence(s) will generate an oligo unless there are too
// few unambiguous bases (in {A,C,G,T}).  Locations that don't
// generate an oligo are silently skipped (except for advancing
// the location counter).
// Auxiliary information:
// -- Whether advancing to this oligo also advanced to a new 
//      sequence within the file. (seqindex == 0)
// -- The contents of the comment/description line for the 
//      current sequence. (descrip)
// -- Whether the normalized oligo returned is the oligo as
//      it occurred in the sequence or its reverse complement.
//      -- istop()
// -- The forward and reverse versions of the oligo, relative
//      to the original sequence.
// -- The location (ending base number) of the oligo in the
//      sequence, with the first base of the sequence numbered 1.
//      -- seqindex
// -- The number of oligos generated in the current sequence,
//      and in the entire file so far.
// -- The number of bases seen in the entire file so far.  (The number
//      seen in the sequence so far, not counting the current
//      oligo, is the current oligo location.)
// -- The number of sequences seen in the entire file so far
//      (including the current sequence).

#ifndef DEFINED_OLIGOSEQ
#define DEFINED_OLIGOSEQ 1
#include "OligoGen.hh"
#include <string.h>
#include <iostream>
#include <fstream>

using namespace std;

class OligoSeq: public OligoGen {
 protected:
  static const int BUFSIZE = 2048;
  char buf[BUFSIZE];
  char descrip[BUFSIZE+1];
  Index bufindex;       // location in buffer
  Index sequences;      // number of descriptions seen so far
  Index64 allbases;     // number of bases in all sequences so far
  Index64 unambiguous;  // number of ACGTacgt in all sequences so far
  Index64 alloligos;    // number of complete oligos seen in all, so far
  Index64 seqindex;     // number of bases in current sequence so far
  istream *in;
  bool softmasked;      // if true, treat lower case as masked

  unsigned char nextBase() {
    // Read character by character if FASTA format, but read in a whole sequence
    // (and it's matching-length quality string) if FASTQ or SCARF.
    unsigned char c;
    while (1) {
      if ((c = buf[bufindex]) < 'A') { // not a base
        if (! c) {
          // end of current buffer contents
          if (in->getline(buf, BUFSIZE)) {
            // successfully got more sequence; keep going
            bufindex = 0;
            continue;
            // if ( (in.rdstate() & ifstream::failbit ) != 0 ) then either EOF or buf filled up
          }
          else {
            // end of file
            buf[0] = '\0';
            // Leave these members unchanged to allow query
            // of last sequence name and length.
            // descrip[0] = 0;
            // seqindex = 0;
            return '\0';
          }
        }
        else if ('>' == c) {
          sequences++;
          bufindex = 0;
          (void) strncpy(descrip, buf, BUFSIZE);
          seqindex = 0;
          buf[0] = '\0';
          return '>';
        }
        else if ('#' == c) {
          // Ignore comment, discarding remaining line contents
          buf[0] = '\0';
          bufindex = 0;
          continue;
        }
        else { // treat as a blank
          bufindex++;
          continue;
        }
      }
      else { // treat as a letter
        allbases++;
        seqindex++;
        bufindex++;
        switch (c) {
        case 'a': case 'c': case 'g': case 't':
          if (softmasked) {
            return '?';
          }
          // else fall through
        case 'A': case 'C': case 'G': case 'T':
          unambiguous++;
          return c;
        default:
          return '?'; // ambiguous base or invalid character
        }
      }
    }
  }

 public:
  inline Index seq_count() { return sequences; }
  inline Index64 base_count() { return allbases; }
  inline Index64 oligo_count() { return alloligos; }
  inline Index64 unambiguous_count() { return unambiguous; }
  inline Index64 get_seqindex() { return seqindex; }
  inline Index64 get_oligostart() { return seqindex + 1 - Length; }
  inline char *get_descrip() { return descrip; }

  OligoSeq(Index tLength, istream &t_in, bool soft = false):
    OligoGen(tLength),
    bufindex(0),
    sequences(0),
    allbases(0),
    alloligos(0),
    unambiguous(0),
    seqindex(0),
    in(&t_in),
    softmasked(soft)
  {
    buf[0] = '\0';
    (void) strncpy(descrip, "NO DESCRIPTION YET", BUFSIZE);
  }
  // Get the next k-mer and return its position in the sequence
  // (1-based sequence index of the last base in the k-mer)
  // (return 0 if we're starting a new sequence)
  int nextPos() {
    unsigned char c;
    while (c = nextBase()) {
      if (c >= 'A') {
        if (advancech(c)) {
          // Successful construction of kmer,
          // return location of *last* base.
          alloligos++;
          return seqindex;
        }
      }
      else if ('?' == c) {
        clear(); // ambiguous character such as X or N, keep going
      }
      else if ('>' == c) {
        clear(); // beginning of a sequence, inform user
        return 0;
      }
      else {
        clear();
        // EOF or ERROR
        break;
      }
    }
    // EOF or ERROR
    clear();
    return -1;
  }
};
#endif
