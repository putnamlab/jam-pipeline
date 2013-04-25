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
// Oligos  class
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
#ifndef DEFINED_OLIGOS
#define DEFINED_OLIGOS 1
#include <cctype>
#include <cstdlib>
#include <iostream>
#include "OligoTools.hh"

class Oligos {
public:
#if (ULONG_MAX <= 0xFFFFFFFFUL)
  typedef unsigned long long Index64; // 64 bits used as an index
  typedef unsigned long long Oligo; // 64 bits
#else
  typedef unsigned long Index64;
  typedef unsigned long Oligo; // 64 bits
#endif
#if (UINT_MAX < 0xFFFFFFFFUL)
  typedef unsigned long Index32;
#else
  typedef unsigned Index32;
#endif
  typedef unsigned long Index;     // 32 or 64 bits, whichever is the native long type
  typedef char OligoString[33];    // big enough for 32 bases + \0

  // Choice of ACGT encoding is so that integer sort on k-mers is
  // the same as alphabetic (lexicographic) sort on k-mers represented 
  // as strings.
  static const Oligo A = (0ULL);
  static const Oligo C = (1ULL);
  static const Oligo G = (2ULL);
  static const Oligo T = (3ULL);
  static const Index BASEBITS = 2;
  static const Oligo BASEMASK = (1ULL << BASEBITS) - 1;

  static const Index OLIGOBITS = sizeof(Oligo)*8;
  static const Index OLIGOBASES = OLIGOBITS/BASEBITS;

  // An ASCII character with A,C,G or T, & 0x6 (to take the second and
  // third least-significant bits) yields A->0, C->1, T->2, G->3
  // (AC*T*G coding).
  //
  // Convert to AC*G*T coding by using the &0x6 bits to shift a
  // byte with binary 10110100 = 0xb4
  static const Oligo BASE_CODING = (0xb4ULL);
  inline Oligo char2base(unsigned char c) {
    return (BASE_CODING >> (c & 0x6)) & BASEMASK;
  }
  static const int BASECHARS = 
  ('a' << (A * 8)) | 
  ('c' << (C * 8)) | 
  ('g' << (G * 8)) | 
  ('t' << (T * 8));
  inline char base2char(Oligo b) {
    return (BASECHARS >> (b * 8)) & 0xFF;
  }

  // More details about Oligo representation
  // Layout within 64 bits:
  // [ unused bits ][ oligo ]
  // 63                     0
  // MSB                  LSB
  // unised bits are 63..(64 - BitsUnused)
  // oligo is in bits (63 - BitsUnused)..0
  const Index Length;
  const Index BitsUnused;
  const Oligo ValMask;     // Used to extract Oligo value

  Oligos(Index tLength):
    Length(tLength),
    BitsUnused(OLIGOBITS - (BASEBITS * Length)),
    ValMask(~0ULL >> BitsUnused)
  { }

  inline Oligo min(Oligo a, Oligo b) { return (a < b)? a : b; }
  inline Oligo max(Oligo a, Oligo b) { return (a > b)? a : b; }

  inline Oligo getOligo(Oligo w) {
    return w & ValMask;
  }
  inline void putOligo(Oligo &w, Oligo oligo) {
    w = (w & ~ValMask) | (oligo & ValMask);
  }

  // Should pass buffer of type OligoString == char[33]
  OligoString ostring;
  inline char *bases(Oligo w) {
    char base = '\0';
    Oligo basebits;
    int i;
    w = getOligo(w);
    // Iterate (Length + 1) times; first time through writes the terminal NULL char
    for (i = Length; i >= 0; i--) {
      ostring[i] = base;
      basebits = w & BASEMASK;
      base = base2char(basebits);
      w >>= BASEBITS;
    }
    return ostring;
  }
  inline char *Bases(Oligo w) {
    char base = '\0';
    Oligo basebits;
    int i;
    w = getOligo(w);
    // Iterate (Length + 1) times; first time through writes the terminal NULL char
    for (i = Length; i >= 0; i--) {
      ostring[i] = base;
      basebits = w & BASEMASK;
      base = toupper(base2char(basebits));
      w >>= BASEBITS;
    }
    return ostring;
  }
  inline int gc(Oligo w) {
    w = getOligo(w);
    Oligo basebits;
    int i, count = 0;
    for (i = Length; i > 0; i--) {
      basebits = w & BASEMASK;
      w >>= BASEBITS;
      if (C == basebits || G == basebits) {
        count++;
      }
    }
    return count;
  }
  inline Oligo OligoCell2RC(Oligo w) {
    return Oligo2RC(getOligo(w));
  }
  inline Oligo Oligo2RC(Oligo w) {
    if (w & ~ValMask) {
      std::cerr << "Dying because of call to Oligo2RC with full Oligo cell, not just kmer, "
                << w << std::endl;
      exit(-1);
    }
    Oligo r = 0;

    for (int i = 0; i < Length; i++) {
      r <<= BASEBITS;
      r |= (~w) & BASEMASK;
      w >>= BASEBITS;
    }
    return r;
  }
  inline Oligo Normalize(Oligo w) {
    Oligo r = Oligo2RC(w);
    return (w < r)? w : r;
  }
  inline Oligo Normalize(Oligo w1, Oligo w2) {
    return (w1 < w2)? w1 : w2;
  }
  inline Oligo mutate(Oligo w,
                      Oligo mask,
                      Index shift) {
    return w ^ (mask << (2*shift));
  }
};
#endif
