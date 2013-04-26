#include "OligoSeq.hh"
#include "OligoHash.hh"
#include "getprime.hh"
#include <string>
#include "gzstream.h"
#include <iomanip>
#include <cctype>
#include <cstdio>

const unsigned NKIDS = 10;

Oligos::Index OptOligoLen;
Oligos::Index OptHashSize;
Oligos::Index OptHashSlicing;
Oligos::Index OptHashSlice;
unsigned OptMin;
unsigned OptMax;
string OptInTable;
bool OptSoftMasking;
bool OptSummary;
bool OptAmbiguous;
string OptPatterns;
string OptPositions;
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
    "   -i {InputTable}  ["<< OptInTable     <<"] Input table: type kmer count bitvector [SNPpos SNPxormask SNPflip kmer count bitvector]\n" <<
    "   -o {OligoLen}    ["<< OptOligoLen    <<"] Number of nucleotide bases per oligo (in otherwords, value of k for k-mer).\n" <<
    "   -H {HashSize}    ["<< OptHashSize    <<"] Number of cells in hash table (can choose based on InputTable).\n" <<
		"   -S {Slicing[:Slice]} [" << OptHashSlicing << ":" << OptHashSlice << "] Slicing factor and slice for unpartnered kmers.\n" <<
		// "   -p {pattern[,pattern]*} [" << OptPatterns << "] Patterns for kmers to be loaded: with '-' for major-minor partners, just one of AA/AP/PA/PP for unpartnered.\n"
		"   -P {position[,position]*} [" << OptPositions << "] SNP positions for major-minor kmers to be loaded (the base position in which the partners differ) relative to first/minor kmer of the pair\n" <<
		"   -a               ["<< OptAmbiguous << "] Turn on loading of kmers with ambiguous SNP partners\n" <<
		"   -s               ["<< OptSummary << "] Turn on printing of summary line - one character per kmer position\n" <<
		"   -x               ["<< OptSoftMasking <<"] Turn on soft masking of input reads (treat lowercase as Ns)\n" <<
    "   -m  {Min}        ["<< OptMin <<"]         Minimum total count in reads for kmers\n" <<
    "   -M  {Max}        ["<< OptMax <<"]         Maximum  ''     ''  ''  ''    '' kmers or partnered kmers\n" <<
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
	// unsigned count: 6;            // Total count of k-mer in all parents & offspring (kept in Info1 of unused bits in hash cell)
	unsigned long inLibs;     // bitvector for presence/absence in libraries
	unsigned unambiguous:  1; // Fields below are meaningful only if this bit is set
	unsigned partnered:    1; // 1 means found, 0 means searched and not found (not unambiguous, as above,
	                          //      can mean too many or too common partners.
	unsigned pos:          5; // offset of SNP base within kmer, 0..31
	unsigned xormask:      2; // mask for SNP base (see Transverse/Transit/Complement above);
	                          //      0=NoMate => no partner kmer found => nonpolymorphic k-locus
	unsigned flip:         1; // set if allelic partner kmer is in table RC relative to this kmer
}; // total bits: 22 (28 with count put back in)

//	'N',          //   Non-polymorphic unpartnered
//	'r',          //   possibly Repetitive, unpaired because of ambiguous partnering
//	'e',          //   -         kmer not found in hash, because not loaded or possibly Error kmer
//	'q',          //   -         no kmer because of Ns or low quality
//	'-',          //   -         no kmer because of slicing
//	'.',          //   -         single . indicates gap between forward and reverse read

static bool usePosition[sizeof(Oligos::Oligo)*4 +1] = { false };
void positionAddList(string posns) {
	// positions comma-separated should be decimal integers in [1..k]
	if (posns[0] == '*') {
		for (int i = 0; i < sizeof(usePosition)/sizeof(bool); i++) {
			usePosition[i] = true;
			if (debugging("o")) cerr << "Added position " << i << " to filter\n";
		}
		return;
	}
	size_t comma;
	for (comma = 0; comma != string::npos; posns = posns.substr(comma + 1)) {
		long num = strtol(posns.c_str(), NULL, 0);
		usePosition[num] = true;
		if (debugging("o")) cerr << "Added position " << num << " from " << posns << " to filter\n";
		comma = posns.find(',');
	}
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

// -p <patterns> 	# patterns: parental presence-absence patterns, comma-separated, to be loaded
// 	e.g. the default of 
// 	PA-PP,AP-PP,AP-PA,AA-PA,AA-PP
// 	(absent in one parent, partner present in both, allowing for kmers having been missed in p7)
// 	(can be improved by number of offspring bits set, especially for uncertain cases like AA-PA)
// 	ALSO includes patterns for non-polymorphic cases, e.g. "PP" for the most common unpartnered case)
// 	NOTE 1: fraction of SNPmer pairs (x for ones that are almost certainly useless, m/p/u for probably useful
// 	x	AA-AA	0.13%
// 		AA-AP	0.56%
// 	u	AA-PA	3.4%
// 	u	AA-PP	6.1%
// 	x	AP-AP	0.48%
// 	u	AP-PA	16.4%
// 	m	AP-PP	19.8%
// 	x	PA-PA	4.0%
// 	p	PA-PP	32.4%
// 	x	PP-PP	16.9%
// 	NOTE 2: patterns are also applied to ambiguous unpartnered kmers, if on input and -ambiguous is set
// 
// -a	# ambiguous: flag to load ambiguously unpartnered kmers
// 
// -mn <number>		# minor allele min total count
// -Mn <number>		# major allele min total count
// -mx <number>		# minor allele max total count
// -Mx <number>		# major allele max total count
// -On <number>		# other kmer min total count
// -Ox <number>		# other kmer max total count
// 	All kmers will be loaded (except as governed by patterns above), but the mins and maxes
// 	(which won't affect the size of table much) will be used in verbose and summary-string reporting
// 
// -s 		# summary: flag; if set print a summary string for each read, one letter per kmer
// 	*	marks the most useful cases for pairs:
// 
// 		a/A	absent: neither allele sampled in either parent (AA-AA)
// 		w/W	only major seen, only in maternal line (w=inverted m) (AA-AP)
// 	*	b/B	only major seen, only in paternal line (b=inverted p) (AA-PA)
// 	*	u/U	unknown-minor-major: unsure of minor-allele parent (AA-PP)
// 		s/S	second: both alleles found only in second parent, p7 (AP-AP)
// 	*	h/H	tied homozygous, ordering arbitrary (AP-PA) (might also be hidden m/M)
// 	*	m/M	maternal-minor-major: minor in p7 (AP-PP)
// 		f/F	first: both alleles found only in first parent, p6 (PA-PA)
// 	*	p/P	paternal-minor-major: minor in p6 (PA-PP)
// 		t/T	tied heterozygous, ordering arbitrary (PP-PP)
// 
// 	q	no kmer due to quality or Ns
//  - for kmers that were sliced out before even checking the table for presence
// 	e	putative error kmer not in table (or omitted for other reasons; see -pa)
// 	[The following other cases might show up as "e" if not loaded into table]
// 	o	kmer's count in reads (or its partner's) was over threshold
// 	r	unpartnered and possibly repetitive; no unique unambiguous kmer partner
// 	N	non-polymorphic unpartnered (reassuringly common case)
// 	.	indicates the gap between forward and reverse read (reverse read
//    reports will be given in order, but position will be on the insert
//    starting with position 300 and counting inwards), and the summary
//    line will be for the whole insert with reverse read flipped as
//    if on the same strand as forward

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
    "                   [$Revision: 1.1 $]\n";
  PrintOptions();
}

//---------------------------------------------------
// * SetupOptions
//---------------------------------------------------
//
int SetupOptions(int argc, char**argv)
{
  // Default values
	OptOligoLen    = 23;        // -o <small_integer>
	OptHashSize    = 0;         // -H <big_integer>
	OptHashSlicing = 11;        // -S <small_prime>[:<hashslice in 0..small_prime-1>]
	OptHashSlice   = 0;         // override with :# on OptHashSlicing
  OptInTable     = "";        // -i <filename>
	//	OptPatterns    = "AA-PA,AA-PP,AP-PA,AP-PP,PA-PP";   // -p <pattern>[,<pattern>]*
	OptPositions   = "3,12,21"; // -P <small_integer>[,<small_integer>]*
	OptAmbiguous   = false;     // -a
	OptMin         = 2;         // -m <num>
	OptMax         = ~0;        // -M <num>
	OptSoftMasking = false;     // -x
	OptSummary     = false;     // -s
  OptDebug = "";              // -d <string>

  // Handle the options...
  int i;
  Oligos::Index prime;
  for (i = 1; i < argc; i++) {
    char prefix = argv[i][0];
    char theOption = argv[i][1];
		char subOption = argv[i][2];
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
      case 'i':
				OptInTable = argv[++i]; break;
			case 'P':
				OptPositions = argv[++i]; 
				break;
			case 'a':
				OptAmbiguous = true; break;
			case 's':
				OptSummary = true; break;
			case 'x':
				OptSoftMasking = true; break;
			case 'm':
				OptMin = strtol(argv[++i], NULL, 0);
				break;
			case 'S':
				parseSlicing(argv[++i]); // sets OptHashSlicing and OptHashSlice
				break;
			case 'M':
				OptMax = strtol(argv[++i], NULL, 0);
				break;
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
	// patternAddList(OptPatterns);
	positionAddList(OptPositions);
  if (debugging("o")) PrintOptions();
  return i;
}

// Return values for readKmerRecord
enum KRecordType {
	UNPAIRED   = '0',
	PAIRED     = '1',
	AMBIGUOUS  = 'x',
	NONMUTUAL  = 'p',
	ENDOFKMERS = 0     // i.e., the null character
};
KRecordType readKmerRecord(istream &in, 
													 Oligos::Oligo &kmer1,
													 Oligos::Index &count1,
													 Oligos::Index &bits1,
													 unsigned &pos,
													 unsigned &xormask,
													 unsigned &flip,
													 Oligos::Oligo &kmer2,
													 Oligos::Index &count2,
													 Oligos::Index &bits2)
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
		char type = *buf;
		if (type == '1') {
			// Partnered kmers: read two kmers with their SNP relationship between
			// Should be minor allele first, major second (check that parent bits are consistent with that)
			// AND if filtering, enforce that the SNP position is middle or INSET bases from kmer end
			parsed = sscanf(buf+1, "%lx %lx %lx %u %u %u %lx %lx %lx", &kmer1, &count1, &bits1, &pos, &xormask, &flip, &kmer2, &count2, &bits2);
			if (parsed != 9) {
				cerr << "readKmerRecord failure on line:\n" << buf;
				exit(-1);
			}
			if (! usePosition[pos]) continue;  // SNPmer pair doesn't have SNP in one of the selected positions
			// if (! patternCheckByBV(bits1, bits2)) continue; // Pair doesn't have right presence in one or both parents
		}
		else if (type == '0') {
			// confirmed unpartnered kmer: save only if parent bits are consistent with requested (default none)
			// AND if filtering, enforce that the kmer is congruent to Slice module SlicingFactor
			parsed = sscanf(buf+1, "%lx %lx %lx", &kmer1, &count1, &bits1);
			if (parsed != 3) {
				cerr << "readKmerRecord failure on line:\n" << buf;
				exit(-1);
			}
			// Only stuff in table if in current slice
			if (OptHashSlicing && (kmer1 % OptHashSlicing != OptHashSlice)) continue;
			// kmer isn't present in right parent or parents
			// if (! patternCheckByBV(bits1)) continue;
		}
		else if (type == '#') {
			// skip comments
			continue;
		}
		else if (type == 'p' || type == 'x') {
			if (! OptAmbiguous) continue;
			// 'x' indicates not uniquely partnered
			// 'p' indicates unique partnering from this kmer unrequited from desired partner

			// ambiguously partnered kmer; save only if specifically asked to and parent bits match
			// patterns requested for "confirmed unpartnered" (case '0' above)
			// AND if filtering, enforce that the kmer is congruent to Slice module SlicingFactor
			parsed = sscanf(buf+1, "%lx %lx %lx", &kmer1, &count1, &bits1);
			if (parsed != 3) {
				cerr << "readKmerRecord failure on line:\n" << buf;
				exit(-1);
			}
			// Only stuff in table if in current slice
			if (OptHashSlicing && (kmer1 % OptHashSlicing != OptHashSlice)) continue;
			// kmer isn't present in right parent or parents
			// if (! patternCheckByBV(bits1)) continue;
		}
		else {
			cerr << "readKmerRecord unrecognized record:\n" << buf;
			exit(-1);
		}
		// We get here if no problems with current kmer line, so return with parsed values
		return (KRecordType) type;
	}
	// Only get here if ran out of lines. Return null character, for no kmer processed.
	return ENDOFKMERS; // NUL character
}

void inline printFields(OligoHash oh,
												Allelic side[],
												OligoHash::Oligo* op) {
	OligoHash::Index oi = op - oh.hash;
	cout 
		<< hex << setw(12) << setfill('0') 
		<< oh.getOligo(*op) 
		<< setw(0) << "\t" 
		<< oh.getInfo1(*op) << "\t" 
		<< side[oi].inLibs << "\t"
		;
}

OligoHash::Index insertOrDie(OligoHash oh,
														 OligoHash::Oligo inmer,
														 OligoHash::Index num) {
	OligoHash::Index index;
	if ((oh.insert(inmer, num) == OligoHash::MISSING) &&
			(oh.lookuploc(inmer, index) == OligoHash::FOUND)) {
		return index;
	}
	else {
		// Should have been missing until we inserted it, then found when looking again!
		cerr << "Failed to insert or find again " << hex << inmer << dec << endl;
		exit(-1);
	}
}

inline OligoSeq::Oligo mutate(OligoSeq::Oligo w,
                              OligoSeq::Oligo mask,
                              unsigned short shift) {
  return w ^ (mask << (2*shift));
}

int main(int argc, char *argv[]) {
	// const OligoSeq::Index p6bit = 2;
	// const OligoSeq::Index p7bit = 1;

  int firstNonOption = SetupOptions(argc, argv);

	// Configure hash table based on input table (on standard input)
	// OptionsFromComments(cin);

	// OptHashSize *= (OptHashSlicing*3)/2;  // Because size was read from one slice, and we plan
	                                      // to read all slices, have to scale up.
	// Note that we create the hash table WITHOUT slicing. The hash table needs to store
	// SNPmers from all slices, so we apply slicing outside of the table to 
	// the other kmers as needed.
  OligoHash oh(OptHashSize, 1, 0, // OptHashSlicing = 1, OptHashSlice = 0 ==> no slicing
							 OptOligoLen);      // Using extra bits only for the count (six bits);
	                                // Also means this code works up to OligoLen=29 without change.
	Allelic *side = (Allelic *) calloc(sizeof(Allelic), OptHashSize);

	// Read in kmers from input table
	OligoSeq::Oligo kmer1, kmer2;
	Oligos::Index total1, bits1, total2, bits2, index1, index2;
	unsigned pos, xormask, flip;
	KRecordType type;
	// const OligoSeq::Oligo KIDBITMASK = ((~0ULL) >> (8*sizeof(Oligos::Index) - NKIDS));
	igzstream inTable(OptInTable.c_str());
	for (type = readKmerRecord(inTable, kmer1, total1, bits1, pos, xormask, flip, kmer2, total2, bits2);
			 (type != ENDOFKMERS);
			 type = readKmerRecord(inTable, kmer1, total1, bits1, pos, xormask, flip, kmer2, total2, bits2)) {

		switch (type) {
			case PAIRED: {
				// Need not check kmer normalization because that rule is universal for OligoHash implementation
				index1 = insertOrDie(oh, kmer1, total1);
				side[index1].inLibs = bits1;
				side[index1].unambiguous = 1;
				side[index1].partnered = 1;
				side[index1].pos = pos;
				side[index1].xormask = xormask;
				side[index1].flip = flip;
				index2 = insertOrDie(oh, kmer2, total2);
				side[index2].inLibs = bits2;
				side[index2].unambiguous = 1;
				side[index2].partnered = 1;
				side[index2].pos = (flip? oh.Length + 1 - pos : pos);
				side[index2].xormask = xormask; // unchanged by flip
				side[index2].flip = flip;
				if (debugging("i") && (kmer1 % 999983)) {
					cerr << "Inserting partners " << oh.Bases(kmer1) << " & " ;
					cerr << oh.Bases(kmer2);
					cerr << dec << " ( " << pos << "," << xormask << (flip? ",-)" : ",+)")<< endl;
				}
			} break;

			case UNPAIRED: {
				index1 = insertOrDie(oh, kmer1, total1);
				side[index1].inLibs = bits1;
				side[index1].unambiguous = 1;
				side[index1].partnered = 0;
			} break;

			case AMBIGUOUS: // fallthrough; treat same as NONMUTUAL
			case NONMUTUAL: {
				index1 = insertOrDie(oh, kmer1, total1);
				side[index1].inLibs = bits1;
				side[index1].unambiguous = 0;
				side[index1].partnered = 0;
			} break;
			default: // everything else should be skipped, until EOF
				cerr << "Unrecognized return case " << type << " from readKmerRecord\n";
				exit(-1);
		}
	}

	// Format of a kmer-report line:
	// type	pos	norm	kmer1	count1	bitvector1	SNPpos	mask	flip	kmer2	count2	bitvector2
	// Columns
	// type   Type of kmer (or pair) according to summary categories above. Lowercase means that kmer in read (first of pair here)
	//        corresponds to minor allele (or at least to the first one in the InputTable pairs, when ordering is arbitrary);
	//        Uppercase means that kmer in read (first of pair here) corresponds to the major allele (or at least second one on input)
	// pos    (decimal) Starting position of kmer1 in read or fragment, with 1 & 101 as beginning & end of forward read;
	//        300 & 200 as beginning and end of reverse read, respectively.
	// normed 0 if normalized kmer1 is the same sense as the read; 1 if normalized is reverse-complement
	// kmer1  Normalized hex-encoded kmer1
	// count1 (hex) Count of kmer1 in whole sequencing project as read from InputTable
	// bitvector1 (hex) Bit vector for kmer1 presence in read libraries (MSB=p6, next p7, next child10, last child1)
	// [remaining fields are only present for awbushmfpt (paired) types, not qeorN]
	// SNPpos index of SNP (one-base difference between kmer1 and kmer2) relative to normalized kmer1 (first base has index=1)
	// mask   xormask for converting normalized kmer1 to possibly-unnormalized version of kmer2
	// flip   0 if kmer2 generated by xormask is already normalized; 1 if normalized kmer2 is reverse-complement relative to kmer1
	// kmer2, count2, bitvector2
	//        as above for kmer1: a normalized kmer2, and its count and bitvector as read from InputTable, in this case relating to 
	//        the major allelic kmer of the pair (or if arbitrary, the second in the table line); 
	//        and SNPpos counts from the right end of normalized kmer2 if flip was set

	int nseqs  = 0;
  int seqset = 1;
  int filearg = 0;

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
      unsigned nonminor = 0, minor = 0;
			string summary = "";
			OligoSeq::Index expectedPos = oh.Length; // i.e., the oligolength
      while ((np = kmers.nextPos()) >= 0) {
        if (np > 0) {
					char typechar = '#';
					for (; expectedPos < np; expectedPos++) {
						summary += 'q';
					}
					expectedPos = np + 1;
          OligoSeq::Oligo w_norm = kmers.current();
					OligoSeq::Index wi;                 // original kmer's index
          if (oh.lookuploc(w_norm, wi) == OligoHash::FOUND) {
						OligoSeq::Index bitvector;
						char fwd = (w_norm == kmers.fwd() ? '0' : '1');
						if (side[wi].partnered) {
							OligoSeq::Index pi;               // kmer partner's index
							OligoSeq::Oligo perturb = mutate(w_norm, side[wi].xormask, oh.Length - side[wi].pos);
							OligoSeq::Oligo partner = oh.Normalize(perturb);
							if (debugging("p") && (w_norm % 999983)) {
								cerr << "Partner for " << oh.Bases(w_norm) << " is " ;
								cerr << oh.Bases(perturb)
										 << " (based on " << dec << " ( " << side[wi].pos << "," << side[wi].xormask << (side[wi].flip? ",-)" : ",+)")
										 << ")";
								if (perturb != partner) {
									cerr << " normalized as " << oh.Bases(partner);
								}
								cerr << endl;
							}
							if (oh.lookuploc(partner, pi) != OligoHash::FOUND) {
								cerr << "Failed to find partner " << hex << partner << " of " << w_norm << endl;
								exit(-1);
							}
							// summary += (typechar = typeCharByBV(side[wi].inLibs, side[pi].inLibs));
							summary += (typechar = '1');
							cout << dec << typechar << "\t" << (np - 22) << "\t" << fwd << "\t";
							printFields(oh, side, oh.hash + wi);
							cout << "\t" << dec << side[wi].pos 
									 << "\t" << side[wi].xormask
									 << "\t" << side[wi].flip
									 << "\t";
							printFields(oh, side, oh.hash + pi);
							cout << endl;
							continue;
						}
						else {
							// not partnered
							OligoSeq::Index count = oh.getInfo1(oh.hash[wi]);
							// OligoSeq::Index pbits = side[wi].inLibs >> NKIDS;
							if (count > OptMax || !side[wi].unambiguous) {
								summary += 'r';
								continue;
							}
							typechar = 'N';

							summary += typechar;
							cout << dec << typechar << "\t" << (np - 22) << "\t" << fwd << "\t";
							printFields(oh, side, oh.hash + wi);
							cout << endl;
							continue;
						}
          }
          else {
            // Not in the minor-or-tied-homozygous allele kmer table
						if (OptHashSlicing && (w_norm % OptHashSlicing != OptHashSlice)) {
							summary += '-';
						}
						else {
							summary += 'e';
						}
					}
        }
        else { // ! np, beginning of a sequence fragment (read or contig, have description line)
					if (summary.length()) {
						cout << "# Summary: " << summary << endl;
						summary = "";
					}
					cout << kmers.get_descrip() << endl;
					expectedPos = oh.Length;
        }
      }
			// Summary for last read
			if (summary.length())
				cout << "# Summary: " << summary << endl;
      // now np < 0
      cerr << "done with " << argv[filearg] << endl;
      inputf.close();
    }
  }

  exit(0);
}
