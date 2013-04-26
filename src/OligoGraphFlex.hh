// OligoGraph class
// Paul Havlak, February 2012
// $Header$
//
#ifndef DEFINED_OLIGOGRAPH
#define DEFINED_OLIGOGRAPH 1
#include <vector>
#include "OligoHash.hh"

typedef enum {STRANDWARDS = 0, INWARDS = 1, OUTWARDS = 2, BACKWARDS = 3} OligoOrient;
// strand == 0 means ">"/top/same-as-normalized-source>, 
//           1 means "<"/bottom/opposite-from-normalized-source
// orientation: source>sink> =strandwards, source><sink =inwards, <source sink> =outwards, <source<sink backwards
//
// orientation = (source_strand << 1) | sink_strand
// strandwards: bin00=0
// inwards:     bin01=1
// outwards:    bin10=2
// backwards:   bin11=3

class OligoEdge {
private:
public: 
	static const unsigned BOGUS = 0;
	unsigned sink;
	unsigned nreads:8;  // number of reads (or templates) supporting the edge
	unsigned orient: 2; // right, out, in, left
	unsigned lib:5;     // library 0 is same-read, otherwise library index for insert size
	unsigned dist:8;    // simple for same-read, start of right-read minus start of left-read
	// Total of 32 + 1+2+5+8+6 = 32+22 bits, will be rounded by compiler up to 8 bytes
	OligoEdge(unsigned Sink) : 
   	sink(Sink)
	{}
	OligoEdge() :
		nreads(0)
	{}
	OligoEdge(unsigned InUse, unsigned Orientation, unsigned Library, unsigned Distance, unsigned NReads) :
		orient(Orientation),
		lib(Library),
		dist(Distance),
		nreads(NReads)
	{}
};
inline OligoOrient make_orient(int leftstrand, int rightstrand) {
	return (OligoOrient) ((leftstrand << 1) | rightstrand);
}

typedef vector<OligoEdge> OligoEdges;

class OligoNode {
	friend class OligoEdge;
public:
	OligoEdges up;   // edges to upstream neighbors, relative to normalized kmer
	OligoEdges down; //   "   to downstream neighbors, "      "     "        "

	static const unsigned NONE = ~0U;

	void add_edge(Oligos::Index sink,
								OligoOrient orient,
								unsigned lib,
								unsigned distance,
								unsigned increment = 1) {
		OligoEdges *edges;
		if (orient == STRANDWARDS || orient == INWARDS) {
			edges = &down; // downstream
		}
		else { // orient in OUTWARDS, BACKWARDS
			edges = &up;   // upstream
		}
		unsigned found = NONE;
		unsigned index;
		if (edges->size() && (*edges)[0].nreads == 0) {
			return; // kmer already had too many distinct edges
		}
		for (index = 0; index < edges->size(); index++) {
			if ((*edges)[index].sink == sink) {
				found = index;
			}
		}
		// Either edges[found] == sink OR
		// append a new edge
		if (found != NONE) {
			(*edges)[found].nreads += increment; // supporting an existing edge
			if (((*edges)[found].orient != orient) ||
					((*edges)[found].dist != distance))
				(*edges)[found].dist = 0; // Marks this edge (and this direction in general) as bogus
  			// Only suitable when walking read-edges, not mate-pair edges
		}
		else { // index == edges->size()
			edges->resize(index + 1);
			(*edges)[index].sink   = sink;
			(*edges)[index].orient = orient;
			(*edges)[index].lib    = lib;
			(*edges)[index].dist   = distance;
			(*edges)[index].nreads = increment;
		}
	}
};

#endif
