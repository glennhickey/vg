#ifndef PILEUP_H
#define PILEUP_H

#include <iostream>
#include <algorithm>
#include <functional>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

// This is a collection of protobuf Pileup records that are indexed
// on their position. Pileups can be merged and streamed, and computed
// from Alignments.  The pileup records themselves are essentially
// protobuf versions of lines in Samtools pileup format.  
class Pileups {
public:
    
    Pileups() {}

    // from protobuf stream
    Pileups(istream& in, bool showp = false);
    
    // copy constructor
    Pileups(const Pileups& other) {
        if (this != &other) {
            for (auto& p : other._pileups) {
                insert(new NodePileup(*p.second));
            }
        }
    }

    // move constructor
    Pileups(Pileups&& other) noexcept {
        _pileups = other._pileups;
        other._pileups.clear();
    }

    // copy assignment operator
    Pileups& operator=(const Pileups& other) {
        Pileups tmp(other);
        *this = move(tmp);
        return *this;
    }

    // move assignment operator
    Pileups& operator=(Pileups&& other) noexcept {
        swap(_pileups, other._pileups);
        other._pileups.clear();
        return *this;
    }

    // delete contents of table
    ~Pileups() {
        clear();
    }
    void clear();

    typedef hash_map<int64_t, NodePileup*> PileupHash;
    
    // This maps from Position to Pileup.
    PileupHash _pileups;

    // write to JSON
    void to_json(ostream& out);
    // read from protobuf
    void load(istream& in);
    // write to protobuf
    void write(ostream& out, uint64_t buffer_size = 1000);

    // apply function to each pileup in table
    void for_each(const function<void(NodePileup&)>& lambda);

    // search hash table for node id
    NodePileup* get(int64_t node_id) {
        auto p = _pileups.find(node_id);
        return p != _pileups.end() ? p->second : NULL;
    }
        
    // get a pileup.  if it's null, create a new one and insert it.
    NodePileup* get_create(int64_t node_id) {
        NodePileup* p = get(node_id);
        if (p == NULL) {
            p = new NodePileup();
            p->set_node_id(node_id);
            _pileups[node_id] = p;
        }
        return p;
    }   

    // insert a pileup into the table. it will be deleted by ~Pileups()!!!
    // return true if new pileup inserted, false if merged into existing one
    bool insert(NodePileup* pileup);

    // create / update all pileups from a single alignment
    void compute_from_alignment(VG& graph, Alignment& alignment);

    // create / update all pileups from an edit (called by above).
    // query stores the current position (and nothing else).  
    void compute_from_edit(NodePileup& pileup, int64_t& node_offset, int64_t& read_offset,
                           const Node& node, const Alignment& alignment,
                           const Mapping& mapping, const Edit& edit);
            
    // move all entries in other object into this one.
    // if two positions collide, they are merged.
    // other will be left empty. this is returned
    Pileups& merge(Pileups& other);

    // merge p2 into p1 and return 1. p2 is left an empty husk
    static BasePileup& merge_base_pileups(BasePileup& p1, BasePileup& p2);

    // merge p2 into p1 and return 1. p2 is lef an empty husk
    static NodePileup& merge_node_pileups(NodePileup& p1, NodePileup& p2);

    // get ith BasePileup record
    static BasePileup* get_base_pileup(NodePileup& np, int64_t offset) {
        assert(offset < np.base_pileup_size());
        return np.mutable_base_pileup(offset);
    }
    static const BasePileup* get_base_pileup(const NodePileup& np, int64_t offset) {
        assert(offset < np.base_pileup_size());
        return &np.base_pileup(offset);
    }

    // get ith BasePileup record, create if doesn't exist
    static BasePileup* get_create_base_pileup(NodePileup& np, int64_t offset) {
        for (int64_t i = np.base_pileup_size(); i <= offset; ++i) {
            np.add_base_pileup();
        }
        return get_base_pileup(np, offset);
    }

    // the bases string in BasePileup doesn't allow random access.  This function
    // will parse out all the offsets of snps, insertions, and deletions
    // into separate arrays, each offset is a pair of indexes in the bases and qualities arrays
    static void parse_base_offsets(const BasePileup& bp,
                                   vector<pair<int, int> >& matchOffsets,
                                   vector<pair<int, int> >& insOffsets,
                                   vector<pair<int, int> >& delOffests);

    // transform case of every character in string
    static void casify(string& seq, bool is_reverse) {
        if (is_reverse) {
            transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
        } else {
            transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        }
    }

    // make the sam pileup style token
    static void make_match(string& seq, int64_t from_length, bool is_reverse) {
        if (seq.length() == 0) {
            seq = string(from_length, is_reverse ? ',' : '.');
        } else {
            casify(seq, is_reverse);
        }
    }
    static void make_insert(string& seq, bool is_reverse) {
        casify(seq, is_reverse);
        stringstream ss;
        ss << "+" << seq.length() << seq; 
        seq = ss.str();
    }
    static void make_delete(string& seq, bool is_reverse) {
        casify(seq, is_reverse);
        stringstream ss;
        ss << "-" << seq.length() << seq;
        seq = ss.str();
    }

    static bool base_equal(char c1, char c2, bool is_reverse) {
        char t1 = ::toupper(c1);
        char t2 = ::toupper(c2);
        return is_reverse ? t1 == reverse_complement(t2) : t1 == t2;
    }

    // get a match base value from a pileup value
    static char extract_match(const BasePileup& bp, int offset) {
        char v = bp.bases()[offset];
        assert(v != '+' && v != '-');
        if (v == ',' || v == '.') {
            return ::toupper(bp.ref_base());
        } else if (::islower(v)) {
            return reverse_complement(::toupper(v));
        }
        return v;
    }
};



}

#endif
