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
// OligoGen class
// Paul Havlak, May 2005
// Based on prior work with k-mer hashes back to summer 2002
// $Header$
//
// The intent is to keep k-mers of length up to 32 bases in 
// "Oligo" data objects of 64 bits, that can be represented 
// in one (possibly doublesized) machine integer on current 
// computers and manipulated primarily by integer arithmetic.
// 
// Mostly manipulate "Oligo" objects that are really just integers.
// Oligos is needed only to carry context and fire off 
// operations that need to know the length of the k-mer.
//
// Here are the mechanisms for building up an Oligo with advancing
// bases.
//
#ifndef DEFINED_OLIGOGEN
#define DEFINED_OLIGOGEN 1
#include "Oligos.hh"
class OligoGen: public Oligos {
public:
  OligoGen(Index tLength):
    Oligos(tLength),
    ngood(0)
  { }

protected:
  Oligo f, r;
  Index ngood;
  unsigned char pred; // last base shifted off, > 3 if previous and/or previous oligos are bad
public:
  inline void clear() {
    ngood = 0;
    // A valid current kmer will have f and r reverse complementary.
    // Zeroing f and r will take k advances to overcome.
    f = r = 0;
  }
  inline bool advancech(char b) {
    return advance(char2base(b));
  }
  inline bool advance(Oligo b) {
    pred = (~r) & BASEMASK; // base being shifted out
    // assert(!(b & ~BASEMASK));
    f <<= BASEBITS; f |= b; f &= ValMask;
    r >>= BASEBITS; r |= (~b) << ((OLIGOBITS - BASEBITS) - BitsUnused); 
    r &= ValMask;
    ngood++;
    return (ngood >= Length);
  }
  inline unsigned char shifted_out() {
    // Gives only the two bits coding for [ACGT]
    return pred;
  }
  inline unsigned char shifted_in() {
    // Gives only the two bits coding for [ACGT]
    return f & BASEMASK;
  }
  inline bool okshiftbyone() {
    return (ngood > Length);
    // i.e., ngood >= Length+1,
    // so that both previous & current oligos had Length good bases
  }

  inline Oligo current() {
    return min(f, r);
  }
  inline Oligo maxcurrent() {
    return max(f, r);
  }
  inline Oligo fwd() {
    return f;
  }
  inline Oligo rev() {
    return r;
  }
  inline bool currOnStrand() {
    return f < r;
  }
};
#endif
