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
// OligoCell objects know how to store extra information in the 
// unused bits of an Oligo integer.
//
// Here are the mechanisms for building up an OligoHash out of 
// OligoCell objects.
//

#ifndef DEFINED_OLIGOHASH
#define DEFINED_OLIGOHASH 1
#include "OligoGen.hh"
#include "OligoCells.hh"
#include <stdlib.h>
#include <iostream>
#include <assert.h>

using namespace std;

template<class T2> class OligoHash: public OligoCells {
public:
  typedef enum { FOUND, MISSING, SLICED, FULL } HashFlag;

  Index64 insertions;
  Index64 distinct;
  Oligo *const hash;
  T2 *const side;
  const Index HashPct;
  const Index Size;
  const Index Slicing;
  const Index Slice;
protected:
  Index primes[13];
  const Index NstepPrimes;
public:
  OligoHash(Index HashSize, Index HashSlicing, Index HashSlice,
            Index OligoLen, Index Info2Len = 0, Index Info3Len = 0) :
    OligoCells(OligoLen, Info2Len, Info3Len),
    Size(HashSize),
    HashPct(HashSize/100),
    Slicing(HashSlicing),
    Slice(HashSlice),
    // Could use C++ new, but not sure if that initializes to zero
    hash((Oligo *) calloc(sizeof(Oligo), HashSize)),
    side((T2 *) calloc(sizeof(T2), HashSize)),
    insertions(0),
    distinct(0),
    NstepPrimes((HashSlicing == (sizeof primes)/sizeof(Index) 
                 ? (HashSlicing - 2)
                 : (sizeof primes)/sizeof(Index)))
    {
      // modulo 3 doesn't work well with powers of 4 (see Knuth)
      assert(NstepPrimes != 3);
      assert(Slicing != 3);
      assert(Slicing > Slice);
      assert(HashSize > 1000);
      static int tprimes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41 };
      for (int i = 0; i < NstepPrimes; i++) {
        if (tprimes[i] == NstepPrimes) {
          primes[i] = 43;
        }
        else if (tprimes[i] == Slicing) {
          primes[i] = 47;
        }
        else {
          primes[i] = tprimes[i];
        }
      }
    }

  inline Oligo* first() {
    Index i;
    for (i = 0; i < Size; i++) {
      if (hash[i]) {
        return &(hash[i]);
      }
    }
    return 0;
  }
  inline Oligo* next(Oligo* current) {
    for (current++; current < hash + Size; current++) {
      if (*current) {
        return current;
      }
    }
    return 0;
  }

protected:
  inline bool inslice(Oligo w) {
    return (w % Slicing == Slice);
  }
public:
  HashFlag lookuploc(Oligo key, Index &location) {
    key = getOligo(key);
    if (! inslice(key)) return SLICED;
    Index start, probe, step, tkey;
    Oligo temp = 0;
    
    // make sure to use hardware integer modulus below
    // if (sizeof(Index) < sizeof(Oligo)) {
    //      tkey  = (Index) (key ^ (key >> 31));
    // }
    // else {
    // tkey = key;
    // }
    probe = start = key % Size;
    if ((temp = hash[probe]) && (getOligo(temp) != key)) {

      // NstepPrimes relatively prime to Slicing
      step = primes[key % NstepPrimes];

      do {
        probe += step;
        probe = (probe >= Size)? (probe -= Size) : probe;
      } while (
               (probe != start) 
               && (temp = hash[probe]) 
               && getOligo(temp) != key);
    }
    location = probe;
    if (! temp) {
      return MISSING;
    }
    else if (getOligo(temp) != key) {
      location = (Index) ~0ULL;
      return FULL;
    }
    else {
      return FOUND;
    }
  }
  HashFlag lookup(Oligo key, Oligo &result) {
    Index loc;
    HashFlag flag = lookuploc(key, loc);

    if (flag == FOUND) {
      result = hash[loc];
      return FOUND;
    }
    else {
      result = 0;
      return flag;
    }
  }
  inline void increment(Oligo &val, Index inc1 = 1, Index inc2 = 0, Index inc3 = 0) {
    putInfo1(val, getInfo1(val) + inc1);
    putInfo2(val, getInfo2(val) + inc2);
    putInfo3(val, getInfo3(val) + inc3);
  }

  HashFlag insert(Oligo key, Index info1inc = 1, Index info2inc = 0, Index info3inc = 0) {
    Index loc;
    Oligo newval;
    key = getOligo(key);
    HashFlag flag = lookuploc(key, loc);

    if (flag == SLICED) {
      return SLICED;
    }
    insertions++;
    if (flag == FULL) {
      OligoString tbases;
      cerr << "\nFull hash on word " << bases(key)
           << ", probe " << loc 
           << ", insertion " << insertions 
           << ", distinct " << distinct << endl;
      // assert(flag != FULL);
      return FULL; // less confusing for optimizer?
    }

    if (flag == MISSING) {
      distinct++;
      if (!(distinct % HashPct)) {
        cerr << "OligoHash is " << dec << distinct / HashPct 
             << " percent full." << endl;
      }
      newval = key;
    }
    else { // flag == FOUND
      newval = hash[loc];
    }
    increment(newval, info1inc, info2inc, info3inc);
    hash[loc] = newval;
    return flag;
  }
  void clear() {
    memset(hash, 0, sizeof(Oligo) * Size);
    insertions = distinct = 0;
  }
  void dump(ostream &os, int minreport) {
    Index i;
    os << " oligolen: " << Length << "\tslicing: " << Slicing << "\tslice: " << Slice << "\tinsertions: " << insertions << "\tdistinct: " << distinct << endl;
    Index countBins[1 << Info1Len][1 << Info2Len];
    memset(countBins, 0, (1 << Info1Len) * (1 << Info2Len) * sizeof(Index));

    Oligo hval;

    for (i = 0; i < Size; i++) {
      hval = hash[i];
      countBins[getInfo1(hval)][getInfo2(hval)]++;
      if (minreport > 0 && getInfo1(hval) >= minreport) {
        OligoString tstring;

        os << bases(getOligo(hval)) << '\t'
           << getInfo1(hval) << '\t'
           << getInfo2(hval) << endl;
      }
    }
    Index j;
    for (i = 0; i < (1 << Info1Len); i++) {
      for (j = 0; j < (1 << Info2Len); j++) {
        if (countBins[i][j]) {
          os << i << '\t' << j << '\t' << countBins[i][j] << endl;
        }
      }
    }
  }
};
#endif
