// Paul Havlak, May 2005
// Based on prior work with k-mer hashes back to summer 2002
// $Header: //CBT-depot/CompBioToolsApps/oligotools/OligoHash.hh#1 $
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

class OligoHash: public OligoCells {
public:
  Index64 insertions;
  Index distinct;
  Oligo *const hash;
  const Index HashPct;
  const Index Size;
  const Index Slicing;
  const Index Slice;
protected:
  Index primes[13];
  const Index NstepPrimes;
public:
  typedef enum { FOUND, MISSING, SLICED, FULL } HashFlag;

  OligoHash(Index HashSize, Index HashSlicing, Index HashSlice,
	    Index OligoLen, Index Info2Len = 0, Index Info3Len = 0) :
    OligoCells(OligoLen, Info2Len, Info3Len),
    Size(HashSize),
    HashPct(HashSize/100),
    Slicing(HashSlicing),
    Slice(HashSlice),
    // Could use C++ new, but not sure if that initializes to zero
    hash((Oligo *) calloc(sizeof(Oligo), HashSize)),
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
    if (sizeof(Index) < sizeof(Oligo)) {
      tkey  = (Index) (key ^ (key >> 31));
    }
    else {
      tkey = key;
    }
    probe = start = tkey % Size;
    if ((temp = hash[probe]) && (getOligo(temp) != key)) {

      // NstepPrimes relatively prime to Slicing
      step = primes[tkey % NstepPrimes];

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
    putInfo1(newval, getInfo1(newval) + info1inc);
    putInfo2(newval, getInfo2(newval) + info2inc);
    putInfo3(newval, getInfo3(newval) + info3inc);
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
