#include "OligoSeq.hh"
#include "OligoHashSide.hh"
#include "getprime.hh"
#include <string>
#include "debugging.hh"
#include "gzstream.h"
#include <iomanip>
#include <cctype>
#include <cstdio>
#include <vector>

const unsigned NKIDS = 10;

// These bit patterns, applied to bases encoded as A=0, C=1, G=2, T=3,
//       map you in the specified way (where "Transverse" means the 
//       transversion that isn't to a complementary base)
const unsigned short NoMate = 0, Transverse = 1, Transit = 2, Complement = 3;
class Allelic
{
public:
	Allelic() : // count(0),
		inLibs(0),
		partnered(0),
		pos(0),
		xormask(NoMate),
		flip(0),
		inFwd(0),
		contig(0),
		upDist(0),
		dnDist(0),
		upFuzzy(0),
		dnFuzzy(0)
  {
  }
	void print(ostream &o) {
		o << dec << "side: " << contig << ":" << contigPos << ":" << contigFlip << " "
			<< hex << inLibs << " " 
			<< partnered << ":" << pos << ":" << xormask << ":" << flip << " "
			<< endl;
	}
	// Mapping to contig
	Oligos::Index32 contig;          // Index of containing kmers contig or scaffold
	Oligos::Index32 contigPos:   30; // Starting base position (1-based) within containing contig or scaffold
	Oligos::Index32 contigFlip:   1; // 0 if same strand sense as contig, 1 if opposite
	// Oligos::Index32 count:				10;     // Total count of k-mer in all parents & offspring
  Oligos::Index32 inLibs: (NKIDS+2); // 00 = neither, 11 = both, etc. (ignore kmers not in any offspring)
  Oligos::Index32 partnered:    1; // 1 means unambiguous partner found; 0 is usually confirmed 
	                                 //    nonpolymorphic, for which next three SNP fields are undefined
  Oligos::Index32 pos:          5; // offset of SNP base within kmer, 0..31
  Oligos::Index32 xormask:      2; // mask for SNP base (see Transverse/Transit/Complement above);
                                   //      0=NoMate => no partner kmer found => nonpolymorphic k-locus
  Oligos::Index32 flip:         1; // set if allelic partner kmer is in table RC relative to this kmer
	Oligos::Index32 inFwd:       1; // Used in read processing only
	// links to upstream and downstream kmers in this contig (note that up & down are relative to kmer, not contig head/tail)
	Oligos::Index32 upDist : 7;       // upstream offset distance to prev kmer (gapsize + kmer length)
	Oligos::Index32 upFlip : 1;      // opposite sense? no=0 or yes=1
	Oligos::Index32 upFuzzy : 1;     // fuzzy because only in terms of templates? no=exact b/c of reads=0, yes=fuzzy=1
	Oligos::Index32 dnDist : 7;     // downstream offset distance to next kmer (gapsize + kmer length)
	Oligos::Index32 dnFlip : 1;    // opposite sense? no=0 or yes=1
	Oligos::Index32 dnFuzzy : 1;   // fuzzy because only in terms of templates? no=exact b/c of reads=0, yes=fuzzy=1
	Oligos::Index32 up;
	Oligos::Index32 down;
}; // total bits: 32 (22, plus 10-bit count)

typedef OligoHash<Allelic> OligoHashPlus;

class Contig {
public:
	Oligos::Index32 type : 8;
	Oligos::Index32 activated : 1;
	Oligos::Index32 flip_first : 1;
	Oligos::Index32 flip_last : 1;
	Oligos::Index32 kfirst;
	Oligos::Index32 klast;
	Contig() :
		activated(0)
	{
	}
	void print(ostream &o) {
		o << "contig: " << ((char) type) << ":" << activated 
			<< endl;
	}
};

class Diagonal {
public:
	Oligos::Index32 contigID;
	short diag; // signed integer diff
	bool anti;
	bool revmate;
	unsigned char rposL;
	unsigned char rposR;
	unsigned char nKmers;
	Diagonal() :
		contigID(0),
		nKmers(0),
		anti(false),
		revmate(false)
	{
	}
	void print(ostream &o) {
		o << "diag: " << contigID << ":" << diag << ":" << anti << revmate
			<< ":" << ((int) rposL)
			<< ":" << ((int) rposR)
			<< ":" << nKmers
			<< endl;
	}
};

string OptKmerContigs;
Oligos::Index32 OptNkmerContigs;
Oligos::Index32 OptOligoLen;
Oligos::Index32 OptInfo1Len;
Oligos::Index OptHashSize;
Oligos::Index32 OptHashSlicing;
Oligos::Index32 OptHashSlice;
bool OptSoftMasking;
string OptDebug;

Debugging debug;

void PrintOptions() {
  cerr << "Option values are:\n"<<
    "   -o {OligoLen}    ["<< OptOligoLen << "] Length of oligos (typically in 16..32).\n" <<
    "   -H {HashSize}    ["<< OptHashSize    <<"] Number of cells in hash table (can choose based on InputTable).\n" <<
		"   -S {Slicing[:Slice]} [" << OptHashSlicing << ":" << OptHashSlice << "] Slicing factor and slice for unpartnered kmers.\n" <<
    "   -x {SoftMasking} ["<< OptSoftMasking <<"] Treat lowercase as masked?\n" <<
		"   -c {KmerContigs} ["<< OptKmerContigs <<"] File with contigs/scaffolds as lists of kmers (or paired SNPmers)\n" <<
		"   -n {nKmerContigs} ["<< OptNkmerContigs << "] Numberof kmer contigs, including singleton kmers not in -c KmerContigs file\n" <<
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
    "GenomeLinkContigs  Read in kmer contigs and kmer lists for reads and construct longer\n" <<
    "                   contigs and scaffolds (currently with no distinction between them)\n" <<
		"                   in the same list-of-kmers format\n" <<
    "                   [Parameter defaults are given in square braces.]\n" <<
    "                   [$Revision: 1.1 $]\n";
  PrintOptions();
}

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
  OptHashSlicing = get_prime(OptHashSlicing);
}

//---------------------------------------------------
// * SetupOptions
//---------------------------------------------------
//
int SetupOptions(int argc, char**argv)
{
  // Default values
	OptHashSize    = 0;         // 'H' <large integer>; will deduce from InputTable if omitted
  OptHashSlicing = 11;        // -S <small_prime>[:<hashslice in 0..small_prime-1>]
  OptHashSlice   = 5;         // override with :# on OptHashSlicing
  OptSoftMasking = false;     // -x
	OptKmerContigs = "";        // -c <string>
	OptNkmerContigs = 0;        // -n <integer>
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
			case 'c':
				OptKmerContigs = argv[++i]; break;
      case 'n':
				OptNkmerContigs = strtoll(argv[++i], NULL, 0); 
				break;
      case 'd':
        OptDebug = argv[++i]; 
				debug.record(OptDebug);
				break;
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
  // if (debug.check('o')) 
		PrintOptions();
  return i;
}

inline Oligos::Index32 insertOrDie(OligoHashPlus &oh,
																	 Oligos::Oligo inmer,
																	 Oligos::Index num) {
  Oligos::Index index;
  if ((oh.insert(inmer, num) == OligoHashPlus::MISSING) &&
      (oh.lookuploc(inmer, index) == OligoHashPlus::FOUND)) {
    return index;
  }
  else {
    // Should have been missing until we inserted it, then found when looking again!
    cerr << "Failed to insert or find again " << hex << inmer << dec << endl;
    exit(-1);
  }
}

static unsigned max_contig = 0;

inline void add_to_contig(OligoHashPlus &oh,
													Contig *contigs,
													Oligos::Index32 c_id,
													Oligos::Index32 k_id)
{
	Allelic *kp = &(oh.side[k_id]);
	kp->contig = c_id;
	max_contig = c_id;
	Contig *cp = &(contigs[c_id]);
	if (! cp->activated) {
		cp->kfirst = cp->klast = k_id;
		cp->flip_first = cp->flip_last = kp->contigFlip;
		cp->activated = 1;
	}
	else {
		// Add to kmers list for this contig
		Oligos::Index32 prev = cp->klast;
		unsigned distance    = kp->contigPos - oh.side[prev].contigPos;
		unsigned diff_strand = (cp->flip_last != kp->contigFlip);
		if (cp->flip_last) {
			oh.side[prev].upDist = distance;
			oh.side[prev].up     = k_id;
			oh.side[prev].upFlip = diff_strand;
		}
		else {
			oh.side[prev].dnDist = distance;
			oh.side[prev].down   = k_id;
			oh.side[prev].dnFlip = diff_strand;
		}
		if (kp->contigFlip) {
			kp->dnDist = distance;
			kp->down   = prev;
			kp->dnFlip = diff_strand;
		}
		else {
			kp->upDist = distance;
			kp->up     = prev;
			kp->upFlip = diff_strand;
		}
		cp->klast     = k_id;
		cp->flip_last = kp->contigFlip;
	}
}

const int BUFSIZE = 2048;
void inputKmerContigs(OligoHashPlus &oh,
											Contig *contigs)
{
  // Read in kmers from input table
  igzstream inKc(OptKmerContigs.c_str());

	OligoSeq::Oligo kmer;
	char buf[BUFSIZE];
	buf[0] = '\0';
	Oligos::Index32 current = 0; // current contig ID, 0 is null contig

	while (inKc.getline(buf, BUFSIZE)) {
		// successfully got a line
		Oligos::Index32 parsed;
		Oligos::Oligo kmer1, kmer2;
		unsigned count1, bits1, count2, bits2, index1, index2;
		unsigned pos, xormask, flip;
		Oligos::Index32 cPosn;   // Contig position of kmer locus
		unsigned cFlip;          // Strand of kmer locus (1 = bottom, opposite from contig)
		char type;
		parsed = sscanf(buf, "%c %u %u %lx %x %x %u %u %u %lx %x %x", 
										&type, // kmer or SNPmer type (ignore, redundant with the bitvector information)
										&cPosn, &cFlip, // kmer mapping to contig
										&kmer1, &count1, &bits1, // kmer, count, & bit vector information for only or representative allele
										&pos, &xormask, &flip,   // position of SNP, xormask for base change, and flipped sense of other allele relative to first
										&kmer2, &count2, &bits2  // 2nd kmer, and its count & bit vector in the reads/libraries
										);
		if (12 == parsed) {
			// we have a SNPmer pair, should be with lesser-encoded kmer first (that we use as representative)
			// Need not check kmer normalization because that rule is universal for OligoHash implementation & saved files
			//   -- (except for saved files tied to reads, in which case kmer in read is first)
			index1 = insertOrDie(oh, kmer1, count1);
			oh.side[index1].inLibs = bits1;
			oh.side[index1].partnered = 1;
			oh.side[index1].pos = pos;
			oh.side[index1].xormask = xormask;
			oh.side[index1].flip = flip;
			oh.side[index1].contigPos  = cPosn;
			oh.side[index1].contigFlip = cFlip;
			// index2 = insertOrDie(oh, kmer2, total2);
			// oh.side[index2].inLibs = bits2;
			// oh.side[index2].unambiguous = 1;
			// oh.side[index2].partnered = 1;
			// oh.side[index2].pos = (flip? oh.Length + 1 - pos : pos);
			// oh.side[index2].xormask = xormask; // unchanged by flip
			// oh.side[index2].flip = flip;
			if (debug.check('i') && !(kmer1 % 999983)) {
				cerr << "Inserting partners " << oh.Bases(kmer1) << " & " ;
				cerr << oh.Bases(kmer2);
				cerr << dec << " (" << pos << "," << xormask << (flip? ",-)" : ",+)")<< endl;
			}
			add_to_contig(oh, contigs, current, index1);
		}
		else if (6 == parsed) {
			// it's an unpartnered kmer (nonpolymorphic)
			index1 = insertOrDie(oh, kmer1, count1);
			oh.side[index1].inLibs = bits1;
			oh.side[index1].partnered = 0;
			oh.side[index1].contigPos  = cPosn;
			oh.side[index1].contigFlip = cFlip;
			if (debug.check('i') && !(kmer1 % 999983)) {
				cerr << "Inserting unpaired " << oh.Bases(kmer1) << endl;
			}
			add_to_contig(oh, contigs, current, index1);
		}
		else {
			// something else, like a contig header
			char cType;
			string cName;
			Oligos::Index32 cID; // contig ID

			parsed = sscanf(buf, ">%c.",
											&cType);

			if (1 == parsed) {
				if ('C' == cType)
					break;
				current++;
				contigs[current].type = cType;
				cerr << dec << "Contig# " << current << buf << endl;
			}
			else if ('#' == *buf) {
				// ignore, a comment
			}
			else {
				cerr << "inputKmerContigs failure on line:\n" << buf;
				exit(-1);
			}
		}
	}
	// Only get here if ran out of lines.
	cerr << "Done with inputKmerContigs\n";
}

inline OligoSeq::Index kidbit(Oligos::Index kid) 
{
	return 1 << (kid - 1);
}

inline Oligos::Oligo mutate(Oligos::Oligo w,
														Oligos::Oligo mask,
														unsigned short shift) {
  return w ^ (mask << (2*shift));
}

// Not universally useful to have a null-index value, but good for here,
// and had better not be 0 (which is valid)
const Oligos::Index NULLINDEX = ~(0UL);

inline Oligos::Index repIndex(OligoHashPlus oh,
															Oligos::Oligo w_norm,
															unsigned &oppstrand) { // can be updated if w_rep kmer not the same as w_norm
  OligoSeq::Index wi;      // to get original kmer's index
  if (oh.lookuploc(w_norm, wi) == OligoHashPlus::FOUND) {
    if (oh.side[wi].partnered) {
      OligoSeq::Index pi;               // kmer partner's index
      OligoSeq::Oligo perturb = mutate(w_norm, oh.side[wi].xormask, oh.Length - oh.side[wi].pos);
      OligoSeq::Oligo partner = oh.Normalize(perturb);
      if (oh.lookuploc(partner, pi) != OligoHashPlus::FOUND) {
        cerr << "Failed to find partner " << hex << partner << " of " << w_norm << endl;
        exit(-1);
      }
      if (partner < w_norm) {
        if (partner != perturb) {
          oppstrand = !oppstrand;
        }
        return pi;
      }
      else
        return wi;
    }
    else {
      // not partnered
      // if (oh.side[wi].unambiguous)
			return wi; // unpartnered => representative index is its own
      // else
			// return NULLINDEX; // don't bother with ambiguously pairable kmers
    }
  }
  else
    return NULLINDEX; // not in the table of interesting kmers
}

Oligos::Index lookupOrAdd(OligoHashPlus &oh,
													Oligos::Oligo kmer,
													Oligos::Index count) {

	Oligos::Index index;
	OligoHashPlus::HashFlag hf = oh.lookuploc(kmer, index);

	if (OligoHashPlus::FOUND == hf) {
		return index;
	}
	else {
		return NULLINDEX;
	}
}

// Return values for readKmerRecord
enum KRecordType {
  UNPAIRED   = '0',
  PAIRED     = '1',
  AMBIGUOUS  = 'x', 
	NONMUTUAL  = 'p',
	ENDOFREAD  = 1,
  ENDOFKMERS = 0     // i.e., the null character
};
inline KRecordType nextReadKmerID(istream &in, 
																	OligoHashPlus &oh,
																	Oligos::Index &k_id,
																	string   &rID,
																	unsigned &rPosn,
																	unsigned &rFlip,
																	bool     &revmate,
																	Oligos::Oligo &snpMer
																	)
// Oligos::Index &count3) 
{
  Oligos::Oligo kmer;
  const int BUFSIZE = 2048;
  char buf[BUFSIZE];
  buf[0] = '\0';

	while (in.getline(buf, BUFSIZE)) {
		// cerr << "INPUT: " << buf;
		// successfully got a line
		Oligos::Index32 parsed;
		Oligos::Oligo kmer1, kmer2;
		unsigned count1, bits1, count2, bits2;
		unsigned pos, xormask, flip;

		char type;
		parsed = sscanf(buf, "%c %u %u %lx %x %x %u %u %u %lx %x %x", 
										&type, // kmer or SNPmer type (ignore, redundant with the bitvector information)
										&rPosn, &rFlip, // kmer mapping to contig
										&kmer1, &count1, &bits1, // kmer, count, & bit vector information for only or representative allele
										&pos, &xormask, &flip,   // position of SNP, xormask for base change, and flipped sense of other allele relative to first
										&kmer2, &count2, &bits2  // 2nd kmer, and its count & bit vector in the reads/libraries
										);
		if (12 == parsed) {
			// We have a SNPmer pair, but will have the one in the read first,
			// so we have to check normalization.
			// Also, only insert if representative is not already in the table (from loading contigs)
			unsigned count, bits;
			snpMer = kmer1;
			if (kmer2 < kmer1) {
				kmer  = kmer2;
				count = count2;
				bits  = bits2;
				rFlip ^= flip;
			}
			else {
				kmer  = kmer1;
				count = count1;
				bits  = bits1;
			}
			// Lookup... later, we'll insert if not already there...
			k_id = lookupOrAdd(oh, kmer, count);
			if (NULLINDEX == k_id)
				continue;
			// But don't do the inserting of kmers for singleton contigs yet...
			if (oh.side[k_id].contig)
				return PAIRED;
			// else add a singleton contig...
			// if not adding a singleton contig, fall through to next iteration
		}
		else if (6 == parsed) {
			// it's an unpartnered kmer (nonpolymorphic)
			k_id = lookupOrAdd(oh, kmer1, count1);
			if (NULLINDEX == k_id)
				continue;
			if (oh.side[k_id].contig) 
				return UNPAIRED;
			// else add a singleton contig... to be implemented
			// if not adding a singleton contig, fall through to next iteration
		}
		else {
			// something else, like a contig header
			char scanbuf[BUFSIZE];
			
			parsed = sscanf(buf, ">%s",
											scanbuf);
			if (1 == parsed) {
				// we're good, fall through to next iteration for kmer
				rID = scanbuf;
				revmate = false;
			}
			else if ('#' == *buf) {
				// ignore, a comment
				revmate = true;
			}
			else {
				cerr << "inputKmerContigs failure on line:\n" << buf;
				exit(-1);
			}
		}
	}
  // Only get here if ran out of lines. Return null character, for no kmer processed.
  return ENDOFKMERS; // NUL character
}

inline void addKmer2Diags(OligoHashPlus &oh,
													vector<Diagonal> &dvec,
													Oligos::Index32 kID,
													unsigned readPos,
													unsigned readFlip,
													bool revmate) {
	Oligos::Index32 cID;
	// cerr << "addKmer2Diags: " << dec << kID << ":" << readPos << ":" << readFlip << " ...";

	if (! (cID = oh.side[kID].contig)) {
		// cerr << "not in a contig" << endl;
		return;
	}
	// xor to classify as same-strand/diagonal or opposite/antidiagonal
	bool  anti = readFlip ^ oh.side[kID].contigFlip;
	short diag = (anti ?
								oh.side[kID].contigPos + readPos :
								// define diagonal this way so that it's more often a positive number
								oh.side[kID].contigPos - readPos);
	unsigned x = dvec.size();
	for (unsigned i = 0; i < x; i++) {
		if (dvec[i].contigID == cID &&
				dvec[i].diag == diag &&
				dvec[i].anti == anti &&
				dvec[i].revmate == revmate) {
			dvec[i].nKmers++;
			dvec[i].rposR = readPos;
			// cerr << "accumulated: " << cID << "/" << diag << "/" << anti << endl;
			return;
		}
	}
	// We got here and need to add a diagonal
	dvec.resize(x + 1);
	dvec[x].contigID = cID;
	dvec[x].rposL    = dvec[x].rposR = readPos;
	dvec[x].diag     = diag;
	dvec[x].anti     = anti;
	dvec[x].revmate  = revmate;
	dvec[x].nKmers   = 1;
	// cerr << "added_new: " << cID << "/" << diag << "/" << anti << endl;
}

int byContig() {
}

int byNkmers() {
}

inline void processKmerDiags(string &prID,
														 OligoHashPlus &oh,
														 vector<Diagonal> &dvec,
														 Contig *contigs,
														 bool matesOL,
														 vector<Oligos::Oligo> &snpF,
														 vector<Oligos::Oligo> &snpR) {
	// cerr << "processKmerDiags: " << dvec.size() << "...";
	unsigned maxFwdPos = 0;
	unsigned maxRevPos = 0;
	unsigned minFwdPos = ~0;
	unsigned minRevPos = ~0;
	for (unsigned i = 0; i < dvec.size(); i++) {
		if (dvec[i].revmate) {
			maxRevPos = (maxRevPos < dvec[i].rposR ? dvec[i].rposR : maxRevPos);
			minRevPos = (minRevPos > dvec[i].rposL ? dvec[i].rposL : minRevPos);
		}
		else {
			maxFwdPos = (maxFwdPos < dvec[i].rposR ? dvec[i].rposR : maxFwdPos);
			minFwdPos = (minFwdPos > dvec[i].rposL ? dvec[i].rposL : minFwdPos);
		}
	}

	// Now we can print read-links and mate-links OR filter them down a bit before 
	// adding them to lists by contig pair linked.
	
	// Output information for read pair, or for just one read if unmated
	// pair: (minFwd,maxFwd;minRev,maxRev)        <- where ';' is replaced by ':' if they share kmers
	// unmated: (minFwd,maxFwd)                   <- even if the read was actually a reverse read
	for (unsigned i = 0; i < dvec.size(); i++) {
		Oligos::Index32 cID = dvec[i].contigID;
		cout << dec
				 << cID << "\t"
				 << prID << dec << '('
				 << (maxFwdPos == 0? 0 : minFwdPos) << ',' << maxFwdPos;
		if (maxRevPos > 0) {
			cout << (matesOL ? ':' : ';')
					 << minRevPos << ',' << maxRevPos;
		}
		cout << ')';
		cout << dvec[i].diag << ":"
				 << (dvec[i].anti ? '-' : '+') 
				 << (dvec[i].revmate ? 'r' : 'f')
				 << '['
				 << ((int) (dvec[i].rposL))
				 << ','
				 << ((int) (dvec[i].rposR))
				 << ']'
				 << '#'
				 << ((int) (dvec[i].nKmers));
		if (snpF.size() + snpR.size()) {
			cout << '(';
			for (unsigned jf = 0; jf < snpF.size(); jf++) {
				cout << hex << snpF[jf] << ',';
			}
			cout << ';';
			for (unsigned jr = 0; jr < snpR.size(); jr++) {
				cout << hex << snpR[jr] << ',';
			}
			cout << ')';
		}
		cout << dec << endl;
	}
}

int main(int argc, char *argv[]) {
	const OligoSeq::Index p6bit = 1 << (NKIDS + 1);
	const OligoSeq::Index p7bit = 1 << NKIDS;

  int firstNonOption = SetupOptions(argc, argv);

  if (1) { // debug.check('s')) {
    cerr << "sizes: Oligo(" << sizeof(Oligos::Oligo) 
         << "), Index(" << sizeof(Oligos::Index)
         << "), Allelic(" << sizeof(Allelic) 
         << ")" << endl << flush;
  }

  OligoHashPlus oh(OptHashSize,
									 1, 0, 
									 // No slicing, that should be handled
									 // by previous processing of input files.
									 OptOligoLen);
	Contig kContigs[OptNkmerContigs];

	inputKmerContigs(oh, kContigs);

	cerr << "pt A\n";

	int nseqs  = 0;
	int seqset = 1;
	int filearg = 0;

	vector<Diagonal> dvec;
	dvec.reserve(200);
	vector<Oligos::Index32> kvec;
	vector<Oligos::Oligo>   snpF;
	vector<Oligos::Oligo>   snpR;

	cerr << "pt B\n";

  for (filearg = firstNonOption; filearg < argc; filearg++) {
    if (!strcmp("/", argv[filearg])) {
      seqset++;
      cerr << "Advancing to read-kmers set #" << seqset << ", starting file: " << argv[filearg+1] << endl;
    }
    else {
      cerr << "Opening read-kmers file " << argv[filearg] << endl;
			igzstream inreads(argv[filearg]);

			string prev_rID = "";
			string rID = "";
			unsigned rPosn = 0;
			unsigned rFlip = 0;
			bool revmate = false;
			Oligos::Index k_id;
			bool matesOL = false;
			KRecordType kType;
			Oligos::Oligo snpMer;

			// Changes here to print all the SNPmers in the read (or read-pair)
			// -- (shared with the contig?) x (fwd or reverse)

			// Assignment here, not comparison; 0 value corresponds to end of kmers
			while (kType = nextReadKmerID(inreads, oh, 
																		k_id,
																		rID, rPosn, rFlip, revmate,
																		snpMer)) {
				if (prev_rID.length() && prev_rID != rID) {
					// Starting a new read, so output hits for old one
					if (dvec.size()) {
						processKmerDiags(prev_rID, oh, dvec, kContigs, matesOL, snpF, snpR);
						// sort diagonal records by contig, retain consistent contigs,
						// add/accumulate summary edges between affected contigs
					}
					dvec.clear();
					// Clear information about fwd-reverse kmer sharing in prep for next read pair
					for (int ki = 0; ki < kvec.size(); ki++) {
						oh.side[kvec[ki]].inFwd = 0;
					}
					kvec.clear();
					matesOL = false;
					snpF.clear();
					snpR.clear();
				}
				// map to contig, increment diagonal
				prev_rID  = rID;
				addKmer2Diags(oh, dvec, k_id, rPosn, rFlip, revmate);
				if (revmate) {
					if (oh.side[k_id].inFwd) {
						matesOL = true;
					}
					if (PAIRED == kType) {
						snpR.push_back(snpMer);
					}
				}
				else {
					oh.side[k_id].inFwd = 1;
					kvec.push_back(k_id);
					if (PAIRED == kType) {
						snpF.push_back(snpMer);
					}
				}
			}
			if (dvec.size()) {
				processKmerDiags(prev_rID, oh, dvec, kContigs, matesOL, snpF, snpR);
				// sort diagonal records by contig, retain consistent contigs,
				// add/accumulate summary edges between affected contigs
			}
			dvec.clear();
			// Clear information about fwd-reverse kmer sharing in prep for next read pair
			for (int ki = 0; ki < kvec.size(); ki++) {
				oh.side[kvec[ki]].inFwd = 0;
			}
			kvec.clear();
			matesOL = false;
			snpF.clear();
			snpR.clear();

      cerr << "done with " << argv[filearg] << endl;
      cout << "# Complete for " << argv[filearg] << endl;
      inreads.close();
    }
  }

  exit(0);
}
