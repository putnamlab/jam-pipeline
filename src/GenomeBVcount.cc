#include "OligoSeq.hh"
#include "OligoHashSide.hh"
#include "getprime.hh"
#include <string>
#include "gzstream.h"
#include <iomanip>
#include <cctype>
#include <cstdio>

Oligos::Index OptOligoLen;
Oligos::Index OptHashSize;
Oligos::Index OptHashSlicing;
Oligos::Index OptHashSlice;
bool OptSoftMasking;
string OptDebug;

bool debugging(const char which[]) {
  if (OptDebug.find('+') != string::npos) return true;

  for (int i = 0; which[i]; i++) {
    if (OptDebug.find(which[i]) != string::npos) return true;
  }
  return false;
}

void PrintOptions() {
  cerr << "Option values are:\n"<<
    "   -o {OligoLen}    ["<< OptOligoLen << "] Length of oligos (typically in 16..32).\n" <<
    "   -H {HashSize}    ["<< OptHashSize    <<"] Number of cells in hash table (can choose based on InputTable).\n" <<
    "   -S {Slicing[:Slice]} [" << OptHashSlicing << ":" << OptHashSlice << "] Slicing factor and slice for unpartnered kmers.\n" <<
    "   -x {SoftMasking} ["<< OptSoftMasking <<"] Treat lowercase as masked?\n" <<
    "   -h               Print help information.\n" <<
    "   -d {DebugString} ["<< OptDebug <<"] Each character indicates a debugging option to turn on;\n" <<
    "                                       '+' indicates turn on all debugging.\n" <<
    "   No more options are allowed after the first argument that can't be parsed as an option.\n" <<
    "   What follows is then a list of file names, possibly separated by a '/' token\n" <<
    "   to delineate boundaries between sequence sets (e.g. genome vs. reads, different\n" <<
    "   read sets, etc.) that will be counted separately.\n" <<
    endl;
}

//---------------------------------------------------
// * PrintHelp
//---------------------------------------------------
//
void PrintHelp()
{
  cerr <<"\n"<<
    "   GenomeBVcount   Reads in sequences and builds k-mers table with counts, bit vectors.\n" <<
    "                   Prints kmers, counts, and bit vectors to STDOUT;\n" <<
    "                          debugging messages and histogram of k-mer frequencies to STDERR.\n" <<
    "                   After options, arguments are FASTA sequence files.\n" <<
    "                          Each set of files separated by blanks and '/' are treated as a sequence set.\n" <<
    "                          For example, if post-option arguments are\n" <<
    "                               set1.fam.gz / set2batch1.fam.gz set2batch2.fam.gz\n" <<
    "                          then the first file's kmers will be counted as being in seqset 1,\n" <<
    "                          the others in seqset 2.\n" <<
    "                   Bits for the lowest-numbered seqsets are rightmost in the presence-absence bitvectors.\n" <<
    "                   [Defaults are given in square braces.]\n" <<
    "                   [$Revision$]\n";
  cerr << OligoToolsCredits;
  PrintOptions();
}

// Convert option value of the form 17:3 into a slicing factor and slice #.
void parseSlicing(char *p) {
  int vals[2] = { 0, 0 };
  int i = 0;
  // Start p one back, so that it can be incremented past delimiter on
  // next loop through
  for (p--; *(++p) && i < 2; i++) {
    int num = strtol(p, &p, 0);
    vals[i] = num;
  }
  // Use larger number for slicing factor, smaller number for slice
  // Note that if only one number is specified, slice defaults to 0
  if (vals[1] > vals[0]) {
    OptHashSlicing = vals[1];
    OptHashSlice = vals[0];
  }
  else {
    OptHashSlicing = vals[0];
    OptHashSlice = vals[1];
  }
  // Double-check that OptHashSlicing is prime
  if (OptHashSlicing > 1) {
    OptHashSlicing = get_prime(OptHashSlicing);
  }
}

//---------------------------------------------------
// * SetupOptions
//---------------------------------------------------
//
int SetupOptions(int argc, char**argv)
{
  // Default values
  OptHashSize    = get_prime(99999);         // -H
  OptHashSlicing = 11;        // -S <small_prime>[:<hashslice in 0..small_prime-1>]
  OptHashSlice   = 0;         // override with :# on OptHashSlicing
  OptSoftMasking = false;     // -x
  OptDebug = "";              // 'd'

  // Handle the options...
  int i;
  Oligos::Index prime;
  for (i = 1; i < argc; i++) {
    char prefix = argv[i][0];
    char theOption = argv[i][1];
    if (prefix == '-') {
      switch(theOption) {
      case 'o':
        OptOligoLen = strtol(argv[++i], NULL, 0);
        break;
      case 'H': {
        OptHashSize = strtoll(argv[++i], NULL, 0); 
        prime = get_prime(OptHashSize);
        if (prime < OptHashSize) {
          cerr << "Hash size " << OptHashSize << " reduced to prime number "
               << prime << "." << endl;
          OptHashSize = prime;
        }
      }
        break;
      case 'S':
        parseSlicing(argv[++i]); // sets OptHashSlicing and OptHashSlice
        break;
      case 'x':
        OptSoftMasking = true; break;
      case 'd':
        OptDebug = argv[++i]; break;
      case 'h':
        PrintHelp(); exit(0); 
        break;
      default: 
        cerr << "Unrecognized option: -" << theOption << "\n";
        PrintHelp(); 
        goto EndOptions;
      }
    }
    else {
      cerr << "End of options, " << (argc - i) << " args left\n";
      break;
    }
  }
 EndOptions:
  if (debugging("o")) PrintOptions();
  return i;
}

inline Oligos::Index kidbit(Oligos::Index kid) 
{
  return 1UL << (kid - 1UL);
}

// An OligoHash table with an extra side array of 64-bit integers that will be used as bit vectors.
typedef OligoHash<Oligos::Index64> OligoHashX;

int main(int argc, char *argv[]) {
  int firstNonOption = SetupOptions(argc, argv);

  if (debugging("b")) {
    cerr << dec << "Sizeof Oligos::Index   = " << sizeof(Oligos::Index) << endl;
    cerr << dec << "Sizeof Oligos::Index64 = " << sizeof(Oligos::Index64) << endl;
    Oligos::Index64 test = 0;
    for (int i = 1; i <= 64; i += 3) {
      test |= kidbit(i);
      cerr << "Setting bit #i: " << dec << i 
           << ", bitvector now = " << hex << test 
           << dec << endl;
    }
  }

  OligoHashX oh(OptHashSize,
                OptHashSlicing, OptHashSlice, 
                OptOligoLen);

  int nseqs  = 0;
  int seqset = 1;
  int filearg = 0;
  long bases = 0;
  long unambiguous = 0;
  long oligos = 0;

  for (filearg = firstNonOption; filearg < argc; filearg++) {
    if (!strcmp("/", argv[filearg])) {
      seqset++;
      cerr << "Advancing to sequence set " << seqset << endl;
    }
    else {
      cerr << "Opening sequence file " << argv[filearg] << endl;
      igzstream inputf(argv[filearg]);

      OligoSeq kmers(OptOligoLen, inputf, OptSoftMasking);
      int np;
      while ((np = kmers.nextPos()) >= 0) {
        if (np > 0) {
          OligoSeq::Oligo w = kmers.current();
          OligoSeq::Index wi;
          OligoHashX::HashFlag hf = oh.lookuploc(w, wi);

          if (hf == OligoHashX::FOUND) {
            oh.side[wi] |= kidbit(seqset);
            oh.increment(oh.hash[wi]);
            oh.insertions++;
          }
          else if (hf == OligoHashX::MISSING) {
            oh.hash[wi] = 0;
            oh.side[wi] = kidbit(seqset);
            oh.putOligo(oh.hash[wi], w);
            oh.increment(oh.hash[wi]);
            oh.insertions++;
            oh.distinct++;
          }
        }
        else { // ! np, end of a sequence fragment (read or contig)
          if (! (++nseqs % 100000)) {
            cerr << "@ " << nseqs << " sequences: " << kmers.get_descrip() << endl;
          }
        }
      }
      // now np < 0
      cerr << "done with " << argv[filearg] << " (np= " << np << " )" << endl;
      if (debugging("s")) {
        cerr << "#" << seqset << "\tbase_count:\t"  << kmers.base_count()        << endl 
             << "#" << seqset << "\tunambiguous:\t" << kmers.unambiguous_count() << endl
             << "#" << seqset << "\toligo_count:\t" << kmers.oligo_count()       << endl
          ;
      }
      bases += kmers.base_count();
      unambiguous += kmers.unambiguous_count();
      oligos += kmers.oligo_count();
      inputf.close();
    }
  }

  // Histogram is count for # of kmers with each frequency.
  // Frequency of each kmer is stored in spare bits of each Oligo object in the hash table (info1).
  // (In fact, nonzero info1 doubles as a sign of non-empty Oligo cell.)
  // Let's max out the histogram at kmer frequency 0x3FFF (16383_10),
  // because we're unlikely to be interested in precise counts higher than that --
  // and we can get them from the kmers detail if needed.
  // (Higher-frequency kmers will be counted as having frequence 0x3FFF.)
  const Oligos::Index MAXFREQp1 = 0x4000UL;
  long histogram[MAXFREQp1] = { 0 };
  Oligos::Oligo* op;

  cout << hex;
  for (op = oh.first(); op; op = oh.next(op)) {
    Oligos::Index index = op - oh.hash;
    Oligos::Index freq = oh.getInfo1(*op);
    if (freq < MAXFREQp1) {
      histogram[freq]++;
    }
    else {
      histogram[MAXFREQp1 - 1]++;
    }
    if (freq < 2)
      continue;
    cout << setw((oh.Length + 1) / 2) << setfill('0') << oh.getOligo(*op) 
         << setw(0) 
         << "\t" << freq
         << "\t" << oh.side[index]
         << endl;
  }
  cout << dec << setw(1) << setfill(' ');
  cerr << "# Histogram:" << dec << endl;
  cerr << "# total_bases:\t"   << bases << endl;
  cerr << "# total_unambig:\t" << unambiguous << endl;
  cerr << "# total_oligos:\t"  << oligos << endl;
  for (long hi = 1; hi < (1 << oh.Info1Len); hi++) {
    if (histogram[hi])
      cerr << "# " << hi << "\t" << histogram[hi] << endl;
  }
  exit(0);
}
