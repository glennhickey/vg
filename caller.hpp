#ifndef CALLER_H
#define CALLER_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "pileup.hpp"

namespace vg {

using namespace std;

// We need to break apart nodes but remember where they came from to update edges
// Wrap all this up in this class.  Store pairs of nodes since we can make up to
// two new nodes for any position (diploid assumption throughout). The first
// pair element is mandatory and second optional, and they both need to be
// same length if they both exist. 
struct NodeDivider {
    // offset in original graph node -> up to 2 nodes in call graph
    typedef map<int, pair<Node*, Node*> > NodeMap;
    // Node id in original graph to map above
    typedef hash_map<int, NodeMap> NodeHash;
    NodeHash index;
    // map given node to offset i of node with id in original graph
    // this function can never handle overlaps (and should only be called before break_end)
    void add_fragment(const Node* orig_node, int offset, Node* subnode);
    // break node if necessary so that we can attach edge at specified side
    // this function wil return NULL if there's no node covering the given location
    pair<Node*, Node*> break_end(const Node* orig_node, VG* graph, int offset, bool left_side);
    // erase everything (but don't free any Node pointers, they belong to the graph)
    void clear();
};

// Super simple variant caller, for now written to get bakeoff evaluation bootstrapped.
// Idea: Idependently process Pileup records, using simple model to make calls that
//       take into account read errors with diploid assumption.  
// Output stored as a graph representing the aligned sample.  
class Caller {
public:

    // log of zero
    static const double Log_zero;
    // heterzygous prior (from r MAQ paper)
    static const double Default_het_prior;
    // minimum size of pileup to call a snp
    static const int Default_min_depth;
    // maximum size of pileup to call a snp
    static const int Default_max_depth;
    // minimum number of reads that support snp required to call it
    static const int Default_min_support;
    // same as min_support, but as fraction of total depth
    static const double Default_min_frac;
    // minimum likelihood to call a snp
    static const double Default_min_likelihood;
    // use this score when pileup is missing quality
    static const char Default_default_quality;
    // use to balance alignments to forward and reverse strand
    static const double Default_max_strand_bias;
    
    Caller(VG* graph,
           double het_prior = Default_het_prior,
           int min_depth = Default_min_depth,
           int max_depth = Default_max_depth,
           int min_support = Default_min_support,
           double min_frac = Default_min_frac,
           double min_likelihood = Default_min_likelihood, 
           bool leave_uncalled = false,
           int default_quality = Default_default_quality,
           double max_strand_bias = Default_max_strand_bias,
           ostream* text_calls = NULL);
    ~Caller();
    void clear();

    // input graph
    VG* _graph;
    // output called graph
    VG _call_graph;
    // optional text file of calls
    ostream* _text_calls;

    // buffer for base calls for each position in the node
    // . = reference
    // - = missing
    typedef pair<string, string> Genotype;
    vector<Genotype> _node_calls;
    vector<double> _node_likelihoods;
    // buffer for current node;
    const Node* _node;
    // max id in call_graph
    int64_t _max_id;
    // link called nodes back to the original graph. needed
    // to figure out where edges go
    NodeDivider _node_divider;
    unordered_set<int64_t> _visited_nodes;
    unordered_set<pair<NodeSide, NodeSide> > _called_edges;

    // used to favour homozygous genotype (r from MAQ paper)
    double _het_log_prior;
    double _hom_log_prior;
    // maximum number of nodes to call before writing out output stream
    int _buffer_size;
    // minimum depth of pileup to call variants on
    int _min_depth;
    // maximum depth of pileup to call variants on
    int _max_depth;
    // min reads supporting snp to call it
    int _min_support;
    // minimum fraction of bases in pileup that nucleotide must have to be snp
    double _min_frac;
    // minimum log likelihood for a snp to be called
    double _min_log_likelihood;
    // don't erase positions that aren't aligned to. so ouput will be original
    // graph plus snps. 
    bool _leave_uncalled;
    // if we don't have a mapping quality for a read position, use this
    char _default_quality;
    // min deviation from .5 in proportion of negative strand reads
    double _max_strand_bias;

    // write the call graph
    void write_call_graph(ostream& out, bool json);

    // call every position in the node pileup
    void call_node_pileup(const NodePileup& pileup);

    // call an edge.  remembering it in a table for the whole graph
    void call_edge_pileup(const EdgePileup& pileup);

    // fill in edges in the call graph (those that are incident to 2 call nodes)
    // and add uncalled nodes (optionally)
    void update_call_graph();
    
    // call position at given base
    void call_base_pileup(const NodePileup& np, int64_t offset);
    
    // Find the top-two bases in a pileup, along with their counts
    void compute_top_frequencies(const BasePileup& bp,
                                 const vector<pair<int, int> >& base_offsets,
                                 string& top_base, int& top_count, int& top_rev_count,
                                 string& second_base, int& second_count, int& second_rev_count);
    
    // compute genotype for a position with maximum prob
    double mp_snp_genotype(const BasePileup& bp,
                           const vector<pair<int, int> >& base_offsets,
                           const string& top_base, const string& second_base,
                           Genotype& mp_genotype);
    // compute likelihood of a genotype samtools-style
    double genotype_log_likelihood(const BasePileup& pb,
                                   const vector<pair<int, int> >& base_offsets,
                                   double g, const string& first, const string& second);

    // write graph structure corresponding to all the calls for the current
    // node.  
    void create_node_calls(const NodePileup& np);

    // create augmented edges (except deletes) connecting the new nodes created by create_node_calls
    void create_augmented_edges(list<Node*> nodes1, list<Node*> nodes2, list<int> offsets1,
                      list<int> offsets2);

    // add edges implied by augmented deletes, breaking nodes they connect if necessary
    void create_delete_edges(list<Node*> nodes1, list<Node*> nodes2, list<int> offsets1,
                             list<int> offsets2);

    // make a path corresponding to a snp in the call grpah
    void create_snp_path(int64_t snp_node, bool secondary_snp);

    // write out vcf-like text calls for a node
    void write_text_calls(const NodePileup& pileup);

    // convert nucleotide base into integer
    static int nidx(char c) {
        c = ::toupper(c);
        switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        }
        return 4;
    }

    // and back
    static char idxn(int i) {
        return "ACGTN"[i];
    }

    // log function that tries to avoid 0s
    static double safe_log(double v) {
        return v == 0. ? Log_zero : ::log(v);
    }

    // tranform ascii phred into probability of error
    static double phred2prob(char q) {
        //q -= 33; dont think we need this after all
        assert(q >= 0);
        // error prob = 10^(-PHRED/10)
        return pow(10., -(double)q / (10.));
    }

    // call missing
    bool missing_call(const Genotype& g) {
        return g.first == "-" &&  g.second == "-";
    }

    // call is reference
    bool ref_call(const Genotype& g) {
        return g.first == "." && (g.second == "." || g.second == "-");
    }

    // call is snp
    bool snp_call(const Genotype& g) {
        return !missing_call(g) && !ref_call(g);
    }

    // classify call as 0: missing 1: reference 2: snp
    int call_cat(const Genotype&g) {
        if (missing_call(g)) {
            return 0;
        } else if (ref_call(g)) {
            return 1;
        }
        return 2;
    }
};



}

#endif
