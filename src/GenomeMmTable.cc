#include "OligoSeq.hh"
#include "OligoHash.hh"
#include "getprime.hh"
#include <string>
#include "gzstream.h"
#include <iomanip>
#include <cctype>
#include <cstdio>

// const unsigned NKIDS = 10;

// These bit patterns, applied to bases encoded as A=0, C=1, G=2, T=3,
//       map you in the specified way (where "Transverse" means the 
//       transversion that isn't to a complementary base)
const unsigned short NoMate = 0, Transverse = 1, Transit = 2, Complement = 3;
class Allelic
{
public:
Allelic() : inLibs(0),
						unambiguous(0),
						partnered(0),
						pos(0),
						xormask(NoMate),
						flip(0)
	{
	}
	// unsigned count: 6;     // Total count of k-mer in all parents & offspring
	unsigned long inLibs;  // used as bitvector
	// unsigned inKids:   NKIDS; // (should have at least one bit set for kids!)
	unsigned unambiguous:  1; // Fields below are meaningful only if this bit is set
	unsigned partnered:    1; // 1 means found, 0 means searched and not found (not unambiguous, as above,
	                          //      can mean too many or too common partners.
	unsigned short pos:          5; // offset of SNP base within kmer, 0..31
	unsigned short xormask:      2; // mask for SNP base (see Transverse/Transit/Complement above);
	                          //      0=NoMate => no partner kmer found => nonpolymorphic k-locus
	unsigned short flip:         1; // set if allelic partner kmer is in table RC relative to this kmer
}; // total bits: 12 (18 with count put back in)

Oligos::Index OptOligoLen;
Oligos::Index OptHashSize;
Oligos::Index OptHashSlicing;
Oligos::Index OptHashSlice;
Oligos::Index OptFilterCount = 254;
string OptInTable;
bool OptSoftMasking;
string OptTag;     // Something to remember this run by
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
    "   -t {Tag}         ["<< OptTag         <<"] Nametag for this analysis run.\n" <<
    "   -f {FilterCount} ["<< OptFilterCount <<"] Ignore input kmers whose total count is greater than this\n" <<
    "   -o {OligoLen}    ["<< OptOligoLen    <<"] Number of nucleotide bases per oligo (in otherwords, value of k for k-mer).\n" <<
    "   -H {HashSize}    ["<< OptHashSize    <<"] Number of cells in hash table (can choose based on InputTable).\n" <<
    "   -h               Print help information.\n" <<
    "   -d {DebugString} ["<< OptDebug <<"] Each character indicates a debugging option to turn on;\n" <<
    "                                       '+' indicates turn on all debugging.\n" <<
    "   [standard input]    Text with input kmers and counts as hex numbers\n" <<
    "                       (comments give oligo length, etc.)\n" <<
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
    "   GenomeMerHist   Reads in sequences and builds k-mers table;\n" <<
    "                   prints out # of distinct mers with each frequency up to a limit.\n"<<
    "                   [Defaults are given in square braces.]\n" <<
    "                   [$Revision: 1.2 $]\n";
  PrintOptions();
}

//---------------------------------------------------
// * SetupOptions
//---------------------------------------------------
//
int SetupOptions(int argc, char**argv)
{
  // Default values
	OptOligoLen    = 23;        // 'o'
	OptHashSize    = 0;         // 'H'; will deduce from InputTable if omitted
  OptFilterCount = 254;       // 'f'
  OptInTable     = "";        // 'i'
  OptTag         = "";        // 't'
  OptDebug = "";              // 'd'

  // Handle the options...
  int i;
  Oligos::Index prime;
  for (i = 1; i < argc; i++) {
    char prefix = argv[i][0];
    char theOption = argv[i][1];
    if (prefix == '-') {
      switch(theOption) {
      case 'f':
				OptFilterCount = strtol(argv[++i], NULL, 0);
				break;
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
      case 'd':
        OptDebug = argv[++i]; break;
      case 'i':
				OptInTable = argv[++i]; break;
      case 't':
				OptTag = argv[++i]; break;
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

void OptionsFromComments(istream &in)
{
	const int BUFSIZE = 2048;
	char buf[BUFSIZE], tag[BUFSIZE];
	unsigned long num;
	buf[0] = '\0';

	if (debugging("o")) {
		cerr << "Parsing additional options from inputTable comments\n";
	}
	while (in.getline(buf, BUFSIZE)) {
		// successfully got a line
		char *bufp = buf;
		while (isspace(*bufp)) bufp++;
	  if (*bufp != '#')
			break;
		int parsed = sscanf(bufp, "%s %lu", tag, &num);

		if (!strcmp(tag, "#oligo_len:"))        OptOligoLen    = num;
		else if (!strcmp(tag, "#slicing_fac:")) OptHashSlicing = num;
		else if (!strcmp(tag, "#slice_#:"))     OptHashSlice   = num;
		else if (!strcmp(tag, "#soft_mask:"))   OptSoftMasking = num;
		else if (!strcmp(tag, "#distinct:")) {
			if (!OptHashSize) {
				// use num distinct as hash size
				OptHashSize = num;
			}
			break; // Handling files from before #end_settings comments were added
		}
		else if (!strcmp(tag, "#end_settings")) break;
		else {} // other comments are skipped, but continue processing
	}
	if (debugging("o")) {
		cerr << "OptOligoLen:    " << OptOligoLen << endl
				 << "OptHashSlicing: " << OptHashSlicing << endl
				 << "OptHashSlice:   " << OptHashSlice << endl
				 << "OptSoftMasking: " << OptSoftMasking << endl
				 << "OptHashSize:    " << OptHashSize << endl
			;
	}
}

inline OligoSeq::Oligo mutate(OligoSeq::Oligo w,
															OligoSeq::Oligo mask,
															unsigned short shift) {
	return w ^ (mask << (2*shift));
}

void findPartner(OligoHash &oh,
								 OligoHash::Oligo w,
								 Allelic &snp)
{
	// Go ahead and find a partner, if possible, regardless of SNP position
	// within kmer. Save complete results, and filter later.
	// Because we're going to have to do most of the work anyway now,
	// to make sure there aren't multiple possible partners (those cases
	// we reject by returning false)
	unsigned i, mask;
	snp.partnered = 0;
	bool found = false;

	for (i = 1; i <= oh.Length; i++) {
		for (mask = Transverse; mask <= Complement; mask++) {
			OligoSeq::Oligo mut = mutate(w, mask, oh.Length - i);
			OligoSeq::Oligo mutnorm = oh.Normalize(mut);
			if (debugging("f")) {
				if (mut != mutnorm) {
					cerr << "Normalized " << oh.Bases(mut);
					cerr << " as " << oh.Bases(mutnorm) << endl;
				}
			}
			OligoSeq::Oligo retval;
			if (oh.lookup(mutnorm, retval) == OligoHash::FOUND) {
				if (!found) {
					found = true;
					snp.partnered = 1;
					snp.pos = i;
					snp.xormask = mask;
					snp.flip = (mut != mutnorm);
					if (debugging("f")) {
						cerr << hex << "Found partner " << mut << "/" << mutnorm
								 << " for " << w << endl;
					}
				}
				else { // Found multiple partners, or a partner that was too common
					snp.unambiguous = 0;
					snp.partnered   = 0;
					return; // unpartnered due to ambiguity ambiguous
				}
			}
		}
	}
	// If *found*, we have a unique mapping from this kmer to a potential partner.
	// ("potential" only because they're not mutually unique)
	// Call that enough confirmation for now.
	// snp.partnered already set...
	snp.unambiguous = 1;  // implies that search was unambiguous: unique partner or none
	return;
}

OligoSeq::Oligo readKmerRecord(istream &in, 
															 Oligos::Index &count1, 
															 Oligos::Index &count2)
// Oligos::Index &count3) 
{
	Oligos::Index parsed;
	OligoSeq::Oligo kmer;
	const int BUFSIZE = 2048;
	char buf[BUFSIZE];
	buf[0] = '\0';

	while (in.getline(buf, BUFSIZE)) {
		// successfully got a line
		char *bufp = buf;
		if (!isxdigit(*buf)) {
			continue; // skip past any lines that don't start with a hex-encoded kmer
		}
		// parsed = sscanf(buf, "%qx %lx %lx %lx", &kmer, &count1, &count2, &count3);
		parsed = sscanf(buf, "%lx %lx %lx", &kmer, &count1, &count2);
		if (parsed != 3) {
			// skip past other problem lines (should probably just die...)
			continue;
		}
		return kmer;
	}
	// Only get here if ran out of lines. Since no kmer should be in file with
	// zero counts, this is legit way to convey end-of-file.
	count1 = count2 = 0;
	return ~0;
}

void inline printFields(OligoHash oh,
												Allelic side[],
												OligoHash::Oligo* op) {
	OligoHash::Index oi = op - oh.hash;
	cout 
		<< hex << setw((OptOligoLen + 1)/2) << setfill('0') 
		<< oh.getOligo(*op) 
		<< setw(0) << "\t" 
		<< oh.getInfo1(*op) << "\t" 
		<< side[oi].inLibs
		;
}

int main(int argc, char *argv[]) {
	const OligoSeq::Index p6bit = 2;
	const OligoSeq::Index p7bit = 1;

  int firstNonOption = SetupOptions(argc, argv);

	// Configure hash table based on input table (on standard input)
	// OptionsFromComments(cin);

	// OptHashSize *= (OptHashSlicing*3)/2;  // Because size was read from one slice, and we plan
	                                      // to read all slices, have to scale up.
  OligoHash oh(OptHashSize, 1, 0, // OptHashSlicing = 1, OptHashSlice = 0 ==> no slicing
							 OptOligoLen);      // Using extra bits only for the count (six bits);
	                                // Also means this code works up to OligoLen=29 without change.
	Allelic *side = (Allelic *) calloc(sizeof(Allelic), OptHashSize);

	// Read in kmers from input table
	OligoSeq::Oligo inmer;
	Oligos::Index bitvector, total, index;
	// const OligoSeq::Oligo KIDBITMASK = ((~0ULL) >> (8*sizeof(Oligos::Index) - NKIDS));
	for (inmer = readKmerRecord(cin, total, bitvector);
			 total;  // zero total from readKmerRecord means EOF
			 inmer = readKmerRecord(cin, total, bitvector)) {

		if (total > 1) {             // ignore the total flukes
			
			if ((oh.insert(inmer, total) == OligoHash::MISSING) &&
					(oh.lookuploc(inmer, index) == OligoHash::FOUND)) {
				side[index].inLibs = bitvector;
			}
			else {
				// Should have been missing until we inserted it, then found when looking again!
				cerr << "Failed to insert or find again " << hex << inmer << dec << endl;
				exit(-1);
			}
			if (debugging("h") && !(inmer % 999983)) {
				cerr << "Inserted " << hex << inmer << " with total " << total 
						 << " and bit vector " << side[index].inLibs
						 << dec << endl;
			}
		}
	}

	// Now iterate throught the hash, finding partners
  Oligos::Oligo* op;
  for (op = oh.first(); op; op = oh.next(op)) {
		Oligos::Index oi = op - oh.hash;
		findPartner(oh, oh.getOligo(*op), side[oi]);
		if (debugging("i") && !(oh.getOligo(*op) % 999983)) {
			cerr << "#Checked partner for " << dec 
					 << oh.Bases(oh.getOligo(*op)) << ": "
					 << side[oi].unambiguous << " "
					 << side[oi].partnered << " "
					 << side[oi].pos << " "
					 << side[oi].xormask << " "
					 << side[oi].flip << "\t";
			// printFields(oh, side, op);
			cerr << endl;
		}
  }
  for (op = oh.first(); op; op = oh.next(op)) {
		// Check and print kmers, if they are mutual unique partners
		// (or just the one kmer, if unambiguous unpartnered;
		//  that is having no non-singleton kmer within one edit)
		OligoSeq::Index oi = (op - oh.hash);
		if (side[oi].unambiguous) {
			if (side[oi].partnered) {
				// Look at oh.getOligo(*op)
				// Construct its partner based on side[op - oh.hash]
				OligoHash::Oligo w, partner;
				w = oh.getOligo(*op);
				partner = mutate(w, side[oi].xormask, oh.Length - side[oi].pos);
				if (side[oi].flip) {
					if (debugging("f")) {
						cerr << "About to normalize " << oh.Bases(partner);
						cerr << ", flip of " << oh.Bases(w) << endl;
					}
					partner = oh.Normalize(partner);
				}
				OligoHash::Index pi;
				if (oh.lookuploc(partner, pi) == OligoHash::FOUND) {
					// Require mutual partnership
					if (side[pi].unambiguous && side[pi].partnered) {
						// Order them:
						if ((side[oi].inLibs  < side[pi].inLibs) // incidentally puts minor before major
								||
								(side[oi].inLibs == side[pi].inLibs && w < partner)) {
							// Print here; otherwise print when we visit the partner
							// partnered==1 kmer1	count1	inLibs	pos	xormask	flip	kmer2	count2	inLibs2
							cout << "1\t"; // partnered
							// This kmer info
							printFields(oh, side, op);
							cout << "\t" << dec
								// pos & xormask of SNP relative to first kmer
									 << side[oi].pos << "\t"
									 << side[oi].xormask << "\t"
									 << side[oi].flip << "\t"
								;
							// Partner kmer info
							printFields(oh, side, oh.hash + pi);
							cout << endl;
						}
					}
					else { // not mutual partners; so print special unpartnered case here "p"
						cout << "p\t"; // partnership not mutual in other direction
						printFields(oh, side, op);
						cout << endl;
					}
				}
				else {
					// error and die, shouldn't happen
					cerr << "Failed to find previously discovered partner " 
							 << oh.Bases(partner) 
							 << " of ";
					cerr << oh.Bases(w)
							 << dec << endl;
					// printFields(oh, side, op);
					cerr << endl;
					exit(-1);
				}
			}
			else { // confirmed no partner
				// partnered=0	kmer	count	inLibs
				cout << "0\t"; // unpartnered unambiguous
				printFields(oh, side, op);
				cout << endl;
			}
		}
		else {
			// ambiguous; print for stats purposes only
			// x=>ambiguous kmer count inLibs
			cout << "x\t"; // ambiguous
			printFields(oh, side, op);
			cout << endl;
		}
  }
  exit(0);
}
