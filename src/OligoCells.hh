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
// OligoCells class
// Paul Havlak, March 2004
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
// OligoCells knows about other data being stored in the unused
// bits of the integer.
//
#ifndef DEFINED_OLIGOCELL
#define DEFINED_OLIGOCELL 1
#include "Oligos.hh"
class OligoCells: public Oligos {
public:
  // Details about OligoCells
  // Layout within 64 bits:
  // [ info1 ][ info2 ][ info2 ][ oligo ]
  // 63                                 0
  // MSB                              LSB
  // info1 is in bits 63..(Info1Shift = 64-Info1Len)
  // info2 is in bits (Info1Shift-1)..(Info2Shift)
  // info3 is in bits (Info2Shift-1)..(Info3Shift = 64-BitsUnused)
  //      (that is, BitsUnused by the oligo)
  // oligo is in bits (63 - BitsUnused)..0
  const Index Info1Len;
  const Index Info1Shift;
  const Oligo Info1Mask;
  const Index Info2Len;
  const Index Info2Shift;
  const Oligo Info2Mask;
  const Index Info3Len;
  const Index Info3Shift;
  const Oligo Info3Mask;

  OligoCells(Index tLength, Index tInfo2Len = 0, Index tInfo3Len = 0): // 1 or 2 added fields
    Oligos(tLength),
    Info1Len(BitsUnused - (tInfo2Len + tInfo3Len)),
    Info1Shift(64 - Info1Len),
    Info1Mask(~0ULL >> (OLIGOBITS - Info1Len)),
    Info2Len(tInfo2Len),
    Info2Shift(64 - (Info1Len + Info2Len)),
    Info2Mask(~0ULL >> (OLIGOBITS - Info2Len)),
    Info3Len(tInfo3Len),
    Info3Shift(64 - BitsUnused),
    Info3Mask(~0ULL >> (OLIGOBITS - Info3Len))
  { }

  inline Index getInfo1(Oligo w) {
    return (w >> Info1Shift) & Info1Mask;
  }
  inline void putInfo1(Oligo &w, Index info) {
    w = (w & ~(Info1Mask << Info1Shift)) | (min(info, Info1Mask) << Info1Shift);
  }
  inline Index getInfo2(Oligo w) {
    return (w >> Info2Shift) & Info2Mask;
  }
  inline void putInfo2(Oligo &w, Index info) {
    w = (w & ~(Info2Mask << Info2Shift)) | (min(info, Info2Mask) << Info2Shift);
  }
  inline Index getInfo3(Oligo w) {
    return (w >> Info3Shift) & Info3Mask;
  }
  inline void putInfo3(Oligo &w, Index info) {
    w = (w & ~(Info3Mask << Info3Shift)) | (min(info, Info3Mask) << Info3Shift);
  }
};
#endif
