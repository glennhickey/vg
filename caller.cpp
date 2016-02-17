#include <cstdlib>
#include <stdexcept>
#include "json2pb.h"
#include "caller.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

const double Caller::Log_zero = (double)-1e100;

// these values pretty arbitrary at this point
const double Caller::Default_het_prior = 0.001; // from MAQ
const int Caller::Default_min_depth = 10;
const int Caller::Default_max_depth = 200;
const int Caller::Default_min_support = 1;
const double Caller::Default_min_frac = 0.25;
const double Caller::Default_min_likelihood = 1e-50;
const char Caller::Default_default_quality = 30;
const double Caller::Default_max_strand_bias = 0.25;

Caller::Caller(VG* graph,
               double het_prior,
               int min_depth,
               int max_depth,
               int min_support,
               double min_frac,
               double min_likelihood, 
               bool leave_uncalled,
               int default_quality,
               double max_strand_bias,
               ostream* text_calls):
    _graph(graph),
    _het_log_prior(safe_log(het_prior)),
    _hom_log_prior(safe_log(.5 * (1. - het_prior))),
    _min_depth(min_depth),
    _max_depth(max_depth),
    _min_support(min_support),
    _min_frac(min_frac),
    _min_log_likelihood(safe_log(min_likelihood)),
    _leave_uncalled(leave_uncalled),
    _default_quality(default_quality),
    _max_strand_bias(max_strand_bias),
    _text_calls(text_calls) {
    _max_id = _graph->max_node_id();
}

// delete contents of table
Caller::~Caller() {
    clear();
}

void Caller::clear() {
    _node_calls.clear();
    _node_likelihoods.clear();
    _call_graph = VG();
    _start_node_map.clear();
    _end_node_map.clear();
    _visited_nodes.clear();
}

void Caller::write_call_graph(ostream& out, bool json) {
    if (json) {
        out << pb2json(_call_graph.graph);
    } else {
        _call_graph.serialize_to_ostream(out);
    }
}

void Caller::call_node_pileup(const NodePileup& pileup) {

    _node = _graph->get_node(pileup.node_id());
    assert(_node != NULL);
    assert(_node->sequence().length() == pileup.base_pileup_size());
    
    _node_calls.clear();
    string def_char = _leave_uncalled ? "." : "-";
    _node_calls.assign(_node->sequence().length(), Genotype(def_char, def_char));
    _node_likelihoods.clear();
    _node_likelihoods.assign(_node->sequence().length(), safe_log(0));

    // todo: parallelize this loop
    // process each base in pileup individually
    for (int i = 0; i < pileup.base_pileup_size(); ++i) {
        if (pileup.base_pileup(i).num_bases() >= _min_depth &&
            pileup.base_pileup(i).num_bases() <= _max_depth) {
            call_base_pileup(pileup, i);
        }
    }

    // add nodes and edges created when making calls to the output graph
    // (_side_map gets updated)
    create_node_calls(pileup);

    // stream out vcf-like text of calls if desired
    if (_text_calls != NULL) {
      write_text_calls(pileup);
    }

    _visited_nodes.insert(_node->id());
}

void Caller::call_edge_pileup(const EdgePileup& pileup) {
    if (pileup.num_reads() >= _min_depth &&
        pileup.num_reads() <= _max_depth) {

        // to do: factor in likelihood?
        Edge edge = pileup.edge();
        _called_edges.insert(NodeSide::pair_from_edge(edge));
    }
}

void Caller::update_call_graph() {
    // if we're leaving uncalled nodes, add'em:
    if (_leave_uncalled) {
        function<void(Node*)> add_node = [&](Node* node) {
            if (_visited_nodes.find(node->id()) == _visited_nodes.end()) {
                _call_graph.create_node(node->sequence(), node->id());
                _start_node_map[node->id()] = NodePair(node->id(), -1);
                _end_node_map[node->id()] = NodePair(node->id(), -1);
            }
        };
        _graph->for_each_node(add_node);
    }

    // map every edge in the original graph to equivalent sides
    // in the call graph. if both sides exist, make an edge in the call graph
    function<void(Edge*)> map_edge = [&](Edge* edge) {
        pair<NodeSide, NodeSide> sides = NodeSide::pair_from_edge(edge);
        // skip uncalled edges if not writing augmented graph
        if (!_leave_uncalled && _called_edges.find(sides) == _called_edges.end()) {
            return;
        }
        const NodeMap& from_map = edge->from_start() ? _start_node_map : _end_node_map;
        auto mappedNode1 = from_map.find(edge->from());
        if (mappedNode1 != from_map.end()) {
            const NodeMap& to_map = edge->to_end() ? _end_node_map : _start_node_map;
            auto mappedNode2 = to_map.find(edge->to());
            if (mappedNode2 != to_map.end()) {
                // up to 2 mappings for each node, so 4 possible edges:
                if (mappedNode1->second.first != -1 && mappedNode2->second.first != -1) {
                    _call_graph.create_edge(mappedNode1->second.first, mappedNode2->second.first,
                                            edge->from_start(), edge->to_end());
                }
                if (mappedNode1->second.first != -1 && mappedNode2->second.second != -1) {
                    _call_graph.create_edge(mappedNode1->second.first, mappedNode2->second.second,
                                            edge->from_start(), edge->to_end());
                }
                if (mappedNode1->second.second != -1 && mappedNode2->second.first != -1) {
                    _call_graph.create_edge(mappedNode1->second.second, mappedNode2->second.first,
                                            edge->from_start(), edge->to_end());
                }
                if (mappedNode1->second.second != -1 && mappedNode2->second.second != -1) {
                    _call_graph.create_edge(mappedNode1->second.second, mappedNode2->second.second,
                                            edge->from_start(), edge->to_end());
                }
            }
        }
    };
    _graph->for_each_edge(map_edge);

    // make sure paths are saved
    _call_graph.paths.rebuild_node_mapping();
    _call_graph.paths.rebuild_mapping_aux();
    _call_graph.paths.to_graph(_call_graph.graph);

    // fix ids
    //_call_graph.sort();
    //_call_graph.compact_ids();
}

void Caller::call_base_pileup(const NodePileup& np, int64_t offset) {
    const BasePileup& bp = np.base_pileup(offset);
    
    // parse the pilueup structure
    vector<pair<int, int> > base_offsets;
    Pileups::parse_base_offsets(bp, base_offsets);

    // compute top two most frequent bases and their counts
    string top_base;
    int top_count;
    int top_rev_count;
    string second_base;
    int second_count;
    int second_rev_count;
    compute_top_frequencies(bp, base_offsets, top_base, top_count, top_rev_count,
                            second_base, second_count, second_rev_count);

    // note first and second base will be upper case too
    string ref_base = string(1, ::toupper(bp.ref_base()));

    // compute threshold
    int min_support = max(int(_min_frac * (double)bp.num_bases()), _min_support);

    // compute strand bias
    double top_sb = top_count > 0 ? abs(0.5 - (double)top_rev_count / (double)top_count) : 0;
    double second_sb = second_count > 0 ? abs(0.5 - (double)second_rev_count / (double)second_count) : 0;
    
    // compute max likelihood snp genotype.  it will be one of the three combinations
    // of the top two bases (we don't care about case here)
    pair<string, string> g;
    _node_likelihoods[offset] = mp_snp_genotype(bp, base_offsets, top_base, second_base, g);

    cerr << "bp " << pb2json(bp) << endl;
    cerr <<"top " << top_base << ":" << top_count << " second " << second_base << ":" << second_count
         << " ml " << _node_likelihoods[offset] << "(minml=" << _min_log_likelihood <<", hp="
         << _het_log_prior << ")" << endl;


    if (_node_likelihoods[offset] >= _min_log_likelihood) {
        cerr << " (1) pass ML ";
        // update the node calls
        if (top_count >= min_support && top_sb <= _max_strand_bias) {
            cerr << " (1) pass thresh ";
            if (g.first != ref_base) {
                _node_calls[offset].first = g.first;
            } else {
                _node_calls[offset].first = ".";
            }
        }
        cerr << g.second << " sc " << second_count  << " ms " << min_support << " sb " << second_sb
             << " maxsb " << _max_strand_bias << endl;
        if (second_count >= min_support && second_sb <= _max_strand_bias) {
            cerr << " (2) pass thresh ";
            if (g.second != ref_base && g.second != g.first) {
                _node_calls[offset].second = g.second;
            } else {
                _node_calls[offset].second = ".";
            }
        }
    }
    cerr << " -- gt = " << _node_calls[offset].first << "," << _node_calls[offset].second << endl;
    if (_node_calls[offset].first == "-" && _node_calls[offset].second != "-") {
        swap(_node_calls[offset].first, _node_calls[offset].second);
    }
}

void Caller::compute_top_frequencies(const BasePileup& bp,
                                     const vector<pair<int, int> >& base_offsets,
                                     string& top_base, int& top_count, int& top_rev_count,
                                     string& second_base, int& second_count, int& second_rev_count) {

    // histogram of pileup entries (base, indel)
    unordered_map<string, int> hist;
    // same thing but just reverse strand (used for strand bias filter)
    unordered_map<string, int> rev_hist;
    
    const string& bases = bp.bases();

    // compute histogram from pileup
    for (auto i : base_offsets) {
        string val = Pileups::extract(bp, i.first);
        if (hist.find(val) == hist.end()) {
            hist[val] = 0;
            rev_hist[val] = 0;
        }
        ++hist[val];

        // val will always be uppcase / forward strand.  we check
        // the original pileup to see if reversed
        if (bases[i.first] == ',' || ::islower(bases[i.first + val.length() - 1])) {
            ++rev_hist[val];
        }
    }

    // find highest occurring string (reference breaks ties)
    top_base.clear();
    top_count = 0;
    string ref_base = string(1, ::toupper(bp.ref_base()));
    if (hist.find(ref_base) != hist.end()) {
        top_base = ref_base;
        top_count = hist[ref_base];
    }   
    for (auto i : hist) {
        if (i.second > top_count) {
            top_base = i.first;
            top_count = i.second;
        }
    }

    // find second highest occurring string (reference breaks ties)
    // todo: do it in same pass as above
    second_base.clear();
    second_count = 0;
    if (top_base != ref_base && hist.find(ref_base) != hist.end()) {
        second_base = ref_base;
        second_count = hist[ref_base];
    }
    for (auto i : hist) {
        if (i.first != top_base && i.second > second_count) {
            second_base = i.first;
            second_count = i.second;
        }
    }
    assert(top_base == "" || top_base != second_base);
    top_rev_count = rev_hist[top_base];
    second_rev_count = rev_hist[second_base];
}

// Estimate the most probable snp genotype
double Caller::mp_snp_genotype(const BasePileup& bp,
                               const vector<pair<int, int> >& base_offsets,
                               const string& top_base, const string& second_base,
                               Genotype& mp_genotype) {

    string ref_base = string(1, ::toupper(bp.ref_base()));
    
    // genotype with 2 reference bases
    double gl = genotype_log_likelihood(bp, base_offsets, 2, top_base, second_base);
    double mp = _hom_log_prior + gl;
    double ml = gl;
    mp_genotype = make_pair(top_base == ref_base ? ref_base : "-",
                            second_base == ref_base ? ref_base : "-");

    // genotype with 1 reference base
    if (top_base != ref_base) {
        double gl = genotype_log_likelihood(bp, base_offsets, 1, top_base, top_base);
        double p = _het_log_prior + gl;
        cerr << " p1 " << p << endl;
        if (p > mp) {
            mp = p;
            ml = gl;
            mp_genotype = make_pair(top_base,
                                    second_base == ref_base ? ref_base : "-");
        }
    }
    if (second_base != ref_base) {
        double gl = genotype_log_likelihood(bp, base_offsets, 1, second_base, second_base);
        double p = _het_log_prior + gl;
        if (p > mp) {
            mp = p;
            ml = gl;
            mp_genotype = make_pair(top_base == ref_base ? ref_base : "-",
                                    second_base);
        }
    }
    // genotype with 0 reference bases
    if (top_base != ref_base && second_base != ref_base) {
        double gl = genotype_log_likelihood(bp, base_offsets, 0, top_base, second_base);
        double p = _het_log_prior + gl;
        cerr << "p0 " << p << endl;
        if (p > mp) {
            mp = p;
            ml = gl;
            mp_genotype = make_pair(top_base, second_base);
        }
    }

    return ml;
}

// This is adapted from Equation 2 (tranformed to log) from
// A statistical framework for SNP calling ... , Heng Li, Bioinformatics, 2011
// http://bioinformatics.oxfordjournals.org/content/27/21/2987.full
double Caller::genotype_log_likelihood(const BasePileup& bp,
                                       const vector<pair<int, int> >& base_offsets,
                                       double g, const string& first, const string& second) {
    // g = number of ref alleles
    // first, second = alt alleles (can be the same if only considering one)
    double m = 2.; // ploidy. always assume two alleles
    double k = (double)base_offsets.size(); // depth

    double log_likelihood = -k * log(m); // 1 / m^k

    const string& bases = bp.bases();
    const string& quals = bp.qualities();
    double perr;
    string ref_base = string(1, ::toupper(bp.ref_base()));

    for (int i = 0; i < base_offsets.size(); ++i) {
        string base = Pileups::extract(bp, base_offsets[i].first);
        char qual = base_offsets[i].second >= 0 ? quals[base_offsets[i].second] : _default_quality;
        perr = phred2prob(qual);
        if (base == ref_base) {
            // ref
            log_likelihood += safe_log((m - g) * perr + g * (1. - perr));
        } else if (base == first || base == second) {
            // alt
            log_likelihood += safe_log((m - g) * (1. - perr) + g * perr);
        }        
        else {
            // just call anything else an error
            log_likelihood += safe_log(2. * perr);
        }
    }

    return log_likelihood;
}


void Caller::create_node_calls(const NodePileup& np) {
    
    int n = _node->sequence().length();
    const string& seq = _node->sequence();
    int cur = 0;
    int cat = call_cat(_node_calls[cur]);

    // we make two tracks of nodes since we'll have up to two
    // parallel calls at any give position (given diploid assumption everywhere)
    // a node at position i in nodes1 will connect to all nodes at i-1 and i+1
    // in both arrays (if they exist)
    list<Node*> nodes1;
    list<Node*> nodes2;
    // reference start position of original node of ith position in above arrays
    // -1 for inserts and gaps
    list<int> offsets1;
    list<int> offsets2;

    // scan calls, merging contiguous reference calls.  only consider
    // ref / snp / inserts on first pass.  
    
    // scan contiguous chunks of a node with same call
    // (note: snps will always be 1-base -- never merged)
    for (int next = 1; next <= n; ++next) {
        int next_cat = next == n ? -1 : call_cat(_node_calls[next]);

        cerr << _node->id() << "." << cur << " " << _node_calls[cur].first << "," <<_node_calls[cur].second  << endl;
        
        // for anything but case where we merge consec. ref/refs
        if (cat == 2 || cat != next_cat) {

            if (cat == 0) {
                cerr << "Missing" << endl;
                // missing? add NULL to both tracks
                nodes1.push_back(NULL);
                offsets1.push_back(cur);
                nodes2.push_back(NULL);
                offsets2.push_back(cur);
            }        
            else if (cat == 1) {
                // add reference
                string new_seq = seq.substr(cur, next - cur);
                Node* node = _call_graph.create_node(new_seq, ++_max_id);
                cerr << "dbl REF " << pb2json(*node) << endl;
                nodes1.push_back(node);
                offsets1.push_back(cur);
                nodes2.push_back(NULL);
                offsets2.push_back(cur);
            }
            else {
                // some mix of reference and alts
                assert(next == cur + 1);
                
                function<void(list<Node*>&, list<Node*>&, list<int>&, list<int>&, string&, string&)>  lambda =
                    [&](list<Node*>& nodes1, list<Node*> nodes2, list<int>& offsets1, list<int>& offsets2,
                        string& call1, string& call2) {
                
                    if (call1 == ".") {
                        // reference base
                        string new_seq = seq.substr(cur, 1);
                        Node* node = _call_graph.create_node(new_seq, ++_max_id);
                        cerr << "REF NODE " << pb2json(*node) << endl;
                        nodes1.push_back(node);
                        offsets1.push_back(cur);
                    }
                    else if (call1[0] != '-' && call1[0] != '+') {
                        // snp base
                        string new_seq = call1;
                        Node* node = _call_graph.create_node(new_seq, ++_max_id);
                        cerr << "SNP NODE (" << call1 << ") " << pb2json(*node) << endl;
                        nodes1.push_back(node);
                        offsets1.push_back(cur);
                    }
                    else if (call1[0] == '+') {
                        // insert
                        int i = 1;
                        for (; call1[i] <= '9' && call1[i] >= '0'; ++i);
                        string new_seq = call1.substr(i, call1.length() - i);
                        Node* node = _call_graph.create_node(new_seq, ++_max_id);
                        cerr << "INS NODE " << pb2json(*node) << endl;
                        nodes1.push_back(node);
                        // offset -1 since insert has no reference coordinate
                        offsets1.push_back(-1);
                        if (call2[0] != '+') {
                            // use null token to make sure we insert relative to 2nd track if
                            // something there
                            nodes2.push_back(NULL);
                            offsets2.push_back(-1);
                        }
                    }
                };

                // apply same logic to both calls, updating opposite arrays
                lambda(nodes1, nodes2, offsets1, offsets2, _node_calls[cur].first, _node_calls[cur].second);
                lambda(nodes2, nodes1, offsets2, offsets1, _node_calls[cur].second, _node_calls[cur].first);
                if (nodes1.size() == nodes2.size() - 1) {
                    nodes1.push_back(NULL);
                    offsets1.push_back(-1);
                }
                if (nodes1.size() == nodes2.size() + 1) {
                    nodes2.push_back(NULL);
                    offsets2.push_back(-1);
                }
            }
            // shift right
            cur = next;
            cat = next_cat;
        }
    }

    create_edges(nodes1, nodes2, offsets1, offsets2);
}
 

void Caller::create_edges(list<Node*> nodes1, list<Node*> nodes2, list<int> offsets1,
                          list<int> offsets2) {
    // add endpoints to map
    NodePair start(nodes1.front() != NULL ? nodes1.front()->id() : -1,
                   nodes2.front() != NULL ? nodes2.front()->id() : -1);
    NodePair end(nodes1.back() != NULL ? nodes1.back()->id() : -1,
                   nodes2.back() != NULL ? nodes2.back()->id() : -1);
    _start_node_map[_node->id()] = start;
    _end_node_map[_node->id()] = end;

        assert(nodes1.size() == nodes2.size());

    list<Node*>::iterator i1 = nodes1.begin();
    list<Node*>::iterator j1 = i1;
    ++j1;
    list<Node*>::iterator i2 = nodes2.begin();
    list<Node*>::iterator j2 = i2;
    ++j2;

    for (int i = 0; i < nodes1.size() - 1; ++i, ++i1, ++i2, ++j1, ++j2) {
        if (*i1 != NULL) {
            if (*j1 != NULL) {
                _call_graph.create_edge(*i1, *j1);
            }
            if (*j2 != NULL) {
                _call_graph.create_edge(*i1, *j2);
            }
        }
        if (*i2 != NULL) {
            if (*j2 != NULL) {
                _call_graph.create_edge(*i2, *j2);
            }
            if (*j1 != NULL) {
                _call_graph.create_edge(*i2, *j1);
            }
        }
    }
}

void Caller::create_snp_path(int64_t snp_node, bool secondary_snp) {

    // for now we don't write secdonary snp, so we have 1 path per *site*
    // and counting paths will give us somethign comparable to snp count
    // from bcftools
    if (!secondary_snp) {
        stringstream name;
        name << "SNP_" << snp_node;

        Mapping mapping;
        Position* pos = mapping.mutable_position();
        // make path that covers node forward with no edits.  not super
        // useful but will use to count snps... 
        pos->set_node_id(snp_node);
        pos->set_offset(0);
        mapping.set_is_reverse(false);
        
        // note: create_path doesn't seem to work.. too rushed to look into
        //list<Mapping>& mappings = _call_graph.paths.create_path(name.str());

        list<Mapping> mappings;
        mappings.push_back(mapping);
        _call_graph.paths._paths.insert(make_pair(name.str(), mappings));
    }
}

void Caller::write_text_calls(const NodePileup& pileup) {
  
  int n = _node->sequence().length();
  assert (_node_calls.size() >= n);
  assert (pileup.node_id() == _node->id());
  assert (pileup.base_pileup_size() >= n);
  const string& seq = _node->sequence();

  const string cat[3] = {"MISSING", "REF", "SNP"};

  for (int i = 0; i < n; ++i) {
      
    // use 1-based coordinates like vcf
    *_text_calls << _node->id() << "\t" << (i + 1) << "\t"
        // reference base
                 << seq[i] << "\t"
        // two comma-separated alternate bases
                 << _node_calls[i].first << "," << _node_calls[i].second << "\t"
        // category
                 << cat[call_cat(_node_calls[i])] << "\t"
        // likelihood
                 << _node_likelihoods[i] << "\t"
        // pileup for debugging
                 << pileup.base_pileup(i).bases() << "\n";
  }  
}

}
