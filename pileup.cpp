#include <cstdlib>
#include <stdexcept>
#include "json2pb.h"
#include "pileup.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

void Pileups::clear() {
    for (auto& p : _node_pileups) {
        delete p.second;
    }
    _node_pileups.clear();

    for (auto& p : _edge_pileups) {
        delete p.second;
    }
    _edge_pileups.clear();
}

void Pileups::to_json(ostream& out) {
    out << "{\"node_pileups\": [";
    for (NodePileupHash::iterator i = _node_pileups.begin(); i != _node_pileups.end();) {
        out << pb2json(*i->second);
        ++i;
        if (i != _node_pileups.end()) {
            out << ",";
        }
    }
    out << "]," << endl << "\"edge_pileups\": [";
    for (EdgePileupHash::iterator i = _edge_pileups.begin(); i != _edge_pileups.end();) {
        out << pb2json(*i->second);
        ++i;
        if (i != _edge_pileups.end()) {
            out << ",";
        }
    }
    out << "]}" << endl;
}

void Pileups::load(istream& in) {
    function<void(Pileup&)> lambda = [this](Pileup& pileup) {
        extend(pileup);
    };
    stream::for_each(in, lambda);
}

void Pileups::write(ostream& out, uint64_t chunk_size) {

    int64_t count = max(_node_pileups.size(), _edge_pileups.size()) / chunk_size;
    if (max(_node_pileups.size(), _edge_pileups.size()) % chunk_size != 0) {
        ++count;
    }

    NodePileupHash::iterator node_it = _node_pileups.begin();
    EdgePileupHash::iterator edge_it = _edge_pileups.begin();
    Pileup pileup;

    // note: this won't work at all in parallel but presumably write
    // is single threaded...
    function<Pileup&(uint64_t)> lambda = [&](uint64_t i) -> Pileup& {
        pileup.clear_node_pileups();
        pileup.clear_edge_pileups();
        for (int j = 0; j < chunk_size && node_it != _node_pileups.end(); ++j, ++node_it) {
            NodePileup* np = pileup.add_node_pileups();
            *np = *node_it->second;
        }
        // unlike for Graph, we don't bother to try to group edges with nodes they attach
        for (int j = 0; j < chunk_size && edge_it != _edge_pileups.end(); ++j, ++edge_it) {
            EdgePileup* ep = pileup.add_edge_pileups();
            *ep = *edge_it->second;
        }

        return pileup;
    };

    stream::write(out, count, lambda);
}

void Pileups::for_each_node_pileup(const function<void(NodePileup&)>& lambda) {
    for (auto& p : _node_pileups) {
        lambda(*p.second);
    }
}

void Pileups::for_each_edge_pileup(const function<void(EdgePileup&)>& lambda) {
    for (auto& p : _edge_pileups) {
        lambda(*p.second);
    }
}

EdgePileup* Pileups::get_edge_pileup(pair<NodeSide, NodeSide> sides) {
    if (sides.first < sides.second) {
        swap(sides.first, sides.second);
    }
    auto p = _edge_pileups.find(sides);
    return p != _edge_pileups.end() ? p->second : NULL;
}
            
// get a pileup.  if it's null, create a new one and insert it.
EdgePileup* Pileups::get_create_edge_pileup(pair<NodeSide, NodeSide> sides) {
    if (sides.first < sides.second) {
        swap(sides.first, sides.second);
    }
    EdgePileup* p = get_edge_pileup(sides);
    if (p == NULL) {
        p = new EdgePileup();
        p->mutable_edge()->set_from(sides.first.node);
        p->mutable_edge()->set_from_start(!sides.first.is_end);
        p->mutable_edge()->set_to(sides.second.node);
        p->mutable_edge()->set_to_end(sides.second.is_end);
        _edge_pileups[sides] = p;
    }
    return p;
}


void Pileups::extend(Pileup& pileup) {
    for (int i = 0; i < pileup.node_pileups_size(); ++i) {
        insert_node_pileup(new NodePileup(pileup.node_pileups(i)));
    }
    for (int i = 0; i < pileup.edge_pileups_size(); ++i) {
        insert_edge_pileup(new EdgePileup(pileup.edge_pileups(i)));
    }
}

bool Pileups::insert_node_pileup(NodePileup* pileup) {
    NodePileup* existing = get_node_pileup(pileup->node_id());
    if (existing != NULL) {
        merge_node_pileups(*existing, *pileup);
        delete pileup;
    } else {
        _node_pileups[pileup->node_id()] = pileup;
    }
    return existing == NULL;
}

bool Pileups::insert_edge_pileup(EdgePileup* pileup) {
    EdgePileup* existing = get_edge_pileup(NodeSide::pair_from_edge(*pileup->mutable_edge()));
    if (existing != NULL) {
        merge_edge_pileups(*existing, *pileup);
        delete pileup;
    } else {
        _edge_pileups[NodeSide::pair_from_edge(*pileup->mutable_edge())] = pileup;
    }
    return existing == NULL;
}

void Pileups::compute_from_alignment(VG& graph, Alignment& alignment) {
    // note to self: alignment.is_reverse() just means the read
    // got reversed when put into sequence.  
    const Path& path = alignment.path();
    int64_t read_offset = 0;
    vector<int> mismatch_counts;
    count_mismatches(graph, path, mismatch_counts);
    // element i = location of rank i in the mapping array
    vector<int> ranks(path.mapping_size() + 1, -1);
    // keep track of read offset of mapping array element i
    vector<int> in_read_offsets(path.mapping_size());
    vector<int> out_read_offsets(path.mapping_size());
    _running_del = NULL;
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        if (graph.has_node(mapping.position().node_id())) {
            const Node* node = graph.get_node(mapping.position().node_id());
            NodePileup* pileup = get_create_node_pileup(node);
            int64_t node_offset = mapping.position().offset();
            in_read_offsets[i] = read_offset;
            for (int j = 0; j < mapping.edit_size(); ++j) {
                const Edit& edit = mapping.edit(j);
                // process all pileups in edit.
                // update the offsets as we go
                compute_from_edit(*pileup, node_offset, read_offset, *node,
                                  alignment, mapping, edit, mismatch_counts);
            }
            out_read_offsets[i] = read_offset - 1;

            // if we're the last base of the read, kill the base pileup
            // if there are too many hanging inserts. 
            if (read_offset == alignment.sequence().length()) {
                int last_offset = node_offset > 0 ? node_offset - 1 : 0;
                filter_end_inserts(*pileup, last_offset, *node);
            }
        }
        int rank = mapping.rank() <= 0 ? i + 1 : mapping.rank();
        if (rank <= 0 || rank >= ranks.size() || ranks[rank] != -1) {
            cerr << "Error determining rank of mapping " << i << " in path " << path.name() << ": "
                 << pb2json(mapping) << endl;
        }
        else {
            ranks[rank] = i;
        }
    }
    cerr << "doing ranks" << endl;
    // loop again over all the edges crossed by the mapping alignment, using
    // the offsets and ranking information we got in the first pass
    for (int i = 2; i < ranks.size(); ++i) {
        int rank1_idx = ranks[i-1];
        int rank2_idx = ranks[i];
        if (rank1_idx > 0 || rank2_idx > 0) {
            auto& m1 = path.mapping(rank1_idx);
            auto& m2 = path.mapping(rank2_idx);
            auto s1 = NodeSide(m1.position().node_id(), (m1.is_reverse() ? false : true));
            auto s2 = NodeSide(m2.position().node_id(), (m2.is_reverse() ? true : false));
            // no quality gives a free pass from quality filter
            char edge_qual = 127;
            if (!alignment.quality().empty()) {
                char from_qual = alignment.quality()[out_read_offsets[rank1_idx]];
                char to_qual = alignment.quality()[in_read_offsets[rank2_idx]];
                edge_qual = min(from_qual, to_qual);
            }
            if (edge_qual >= _min_quality) {
                EdgePileup* edge_pileup = get_create_edge_pileup(pair<NodeSide, NodeSide>(s1, s2));
                cerr << "adding ep " << s1 << " --> " << s2 << endl;
                edge_pileup->set_num_reads(edge_pileup->num_reads() + 1);
                if (!alignment.quality().empty()) {
                    *edge_pileup->mutable_qualities() += edge_qual;
                }
            }
        }
    }
    
    assert(alignment.sequence().empty() ||
           alignment.path().mapping_size() == 0 ||
           read_offset == alignment.sequence().length());

}

void Pileups::compute_from_edit(NodePileup& pileup, int64_t& node_offset,
                                int64_t& read_offset,
                                const Node& node, const Alignment& alignment,
                                const Mapping& mapping, const Edit& edit,
                                const vector<int>& mismatch_counts) {
    string seq = edit.sequence();
    // is the mapping reversed wrt read sequence? use for iterating
    bool map_reverse = mapping.is_reverse();
    // is the mapping reversed wrt to the *graph*? use for flipping stuff when writing pileup
    bool aln_reverse = alignment.is_reverse() != mapping.is_reverse();
    // is seq reversed compared to what we actually want to write
    bool seq_reverse = alignment.is_reverse();

    // reverse mapping, we flip sequence for comparison purposes to read
    if (map_reverse) {
        seq = reverse_complement(seq);
    }
    
    // ***** MATCH *****
    if (edit.from_length() == edit.to_length()) {
        _running_del = NULL;
        assert (edit.from_length() > 0);
        make_match(seq, edit.from_length(), aln_reverse);
        assert(seq.length() == edit.from_length());            
        int64_t delta = map_reverse ? -1 : 1;
        for (int64_t i = 0; i < edit.from_length(); ++i) {
            if (pass_filter(alignment, read_offset, mismatch_counts)) {
                BasePileup* base_pileup = get_create_base_pileup(pileup, node_offset);
                // reference_base if empty
                if (base_pileup->num_bases() == 0) {
                    base_pileup->set_ref_base(node.sequence()[node_offset]);
                } else {
                    assert(base_pileup->ref_base() == node.sequence()[node_offset]);
                }
                // add base to bases field (converting to ,. if match)
                char base = seq[i];
                if (base_equal(seq[i], node.sequence()[node_offset], false)) {
                    base = aln_reverse ? ',' : '.';
                } else if (seq_reverse && base != ',' && base != '.') {
                    base = reverse_complement(::toupper(base));
                    if (aln_reverse) {
                        base = ::tolower(base);
                    }
                }
                *base_pileup->mutable_bases() += base;
                // add quality if there
                if (!alignment.quality().empty()) {
                    *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
                }
                // pileup size increases by 1
                base_pileup->set_num_bases(base_pileup->num_bases() + 1);
            }
            // move right along read, and left/right depending on strand on reference
            node_offset += delta;
            ++read_offset;
        }
    }
    // ***** INSERT *****
    else if (edit.from_length() < edit.to_length()) {
        _running_del = NULL;
        if (pass_filter(alignment, read_offset, mismatch_counts)) {
            if (seq_reverse) {
                reverse_complement(seq);
            }
            make_insert(seq, aln_reverse);
            assert(edit.from_length() == 0);
            // we define insert (like sam) as insertion between current and next
            // position (on forward node coordinates). this means an insertion before
            // offset 0 is invalid! 
            int64_t insert_offset =  map_reverse ? node_offset : node_offset - 1;
            if (insert_offset >= 0) {        
                BasePileup* base_pileup = get_create_base_pileup(pileup, insert_offset);
                // reference_base if empty
                if (base_pileup->num_bases() == 0) {
                    base_pileup->set_ref_base(node.sequence()[insert_offset]);
                } else {
                    assert(base_pileup->ref_base() == node.sequence()[insert_offset]);
                }
                // add insertion string to bases field
                base_pileup->mutable_bases()->append(seq);
                if (!alignment.quality().empty()) {
                    *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
                }
                // pileup size increases by 1
                base_pileup->set_num_bases(base_pileup->num_bases() + 1);
            }
            else {
                // need to check with aligner to make sure this doesn't happen, ie
                // inserts would hang off the end of previous node instead of start
                // of this node
                /*
                  stringstream ss;
                  ss << "Warning: pileup does not support insertions before 0th base in node."
                  << " Offending edit: " << pb2json(edit) << endl;
                  #pragma omp critical(cerr)
                  cerr << ss.str();
                */
            }
        }
        // move right along read (and stay put on reference)
        read_offset += edit.to_length();
    }
    // ***** DELETE *****
    else {
        if (pass_filter(alignment, read_offset, mismatch_counts)) {
            assert(edit.to_length() == 0);
            assert(edit.sequence().empty());
            int64_t del_start = !map_reverse ? node_offset :
                node_offset - edit.from_length() + 1;
            seq = node.sequence().substr(del_start, edit.from_length());
            // add deletion string to bases field
            if (seq_reverse) {
                reverse_complement(seq);
            }            
            make_delete(seq, node.id(), node_offset, map_reverse);
            if (_running_del != NULL) {
                // we are appending onto existing entry
                // deletes are special in that they can span multiple nodes/edits
                append_delete(*_running_del->mutable_bases(), seq);
            } else {
                BasePileup* base_pileup = get_create_base_pileup(pileup, node_offset);
                _running_del = base_pileup;

                // reference_base if empty
                if (base_pileup->num_bases() == 0) {
                    base_pileup->set_ref_base(node.sequence()[node_offset]);
                } else {
                    assert(base_pileup->ref_base() == node.sequence()[node_offset]);
                }
                base_pileup->mutable_bases()->append(seq);
                if (!alignment.quality().empty()) {
                    *base_pileup->mutable_qualities() += alignment.quality()[read_offset];
                }
                // pileup size increases by 1
                base_pileup->set_num_bases(base_pileup->num_bases() + 1);                
            }
        }
        int64_t delta = map_reverse ? -edit.from_length() : edit.from_length();
        // stay put on read, move left/right depending on strand on reference
        node_offset += delta;
    }
}

void Pileups::count_mismatches(VG& graph, const Path& path,
                               vector<int>& mismatches,
                               bool skipIndels)
{
    mismatches.clear();
    int64_t read_offset = 0;
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        if (graph.has_node(mapping.position().node_id())) {
            const Node* node = graph.get_node(mapping.position().node_id());
            int64_t node_offset = mapping.position().offset();
            for (int j = 0; j < mapping.edit_size(); ++j) {
                const Edit& edit = mapping.edit(j);
                // process all pileups in edit.
                // update the offsets as we go
                string seq = edit.sequence();
                bool is_reverse = mapping.is_reverse();
    
                // ***** MATCH *****
                if (edit.from_length() == edit.to_length()) {
                    int64_t delta = is_reverse ? -1 : 1;
                    for (int64_t i = 0; i < edit.from_length(); ++i) {
                        if (!edit.sequence().empty() &&
                            !base_equal(seq[i], node->sequence()[node_offset], is_reverse)) {
                            mismatches.push_back(1);
                        }
                        else {
                            mismatches.push_back(0);
                        }
                        // move right along read, and left/right depending on strand on reference
                        node_offset += delta;
                        ++read_offset;
                    }
                }
                // ***** INSERT *****
                else if (edit.from_length() < edit.to_length()) {
                    if (skipIndels == false) {
                        mismatches.push_back(1);
                        for (int x = 1; x < edit.to_length(); ++x) {
                            mismatches.push_back(0);
                        }
                    }
                    // move right along read (and stay put on reference)
                    read_offset += edit.to_length();
                }
                // ***** DELETE *****
                else {
                    if (skipIndels == false) {
                        // since we're working in read coordinates, we count
                        // a single mismatch right before the delete.
                        if (mismatches.size() > 0) {
                            mismatches[mismatches.size() - 1] = 1;
                        }
                    }
                    int64_t delta = is_reverse ? -edit.from_length() : edit.from_length();
                    // stay put on read, move left/right depending on strand on reference
                    node_offset += delta;
                }
            }
        }
    }
    assert(skipIndels || read_offset == mismatches.size());
    // too lazy to do full count inline.  sum up here
    int count = 0;
    for (int i = 0; i < mismatches.size(); ++i) {
        count += mismatches[i];
        mismatches[i] = count;
    }
}

bool Pileups::pass_filter(const Alignment& alignment, int read_offset,
                          const vector<int>& mismatches) const
{
    bool passes = true;
    if (!alignment.quality().empty()) {
        passes = alignment.quality()[read_offset] >= _min_quality;
    }
    if (_window_size > 0 && passes) {
        // counts in left window
        int left_point = max(0, read_offset - _window_size / 2 - 1);
        int right_point = max(0, read_offset - 1);
        int count = mismatches[right_point] - mismatches[left_point];
        // coutns in right window
        left_point = read_offset;
        right_point = min(read_offset + _window_size / 2, (int)mismatches.size() - 1);
        count += mismatches[right_point] - mismatches[left_point];
        passes = passes && count <= _max_mismatches;
    }
    return passes;
}

void Pileups::filter_end_inserts(NodePileup& pileup, int64_t node_offset, const Node& node)
{
    int insert_count = 0;
    BasePileup* base_pileup = get_create_base_pileup(pileup, node_offset);
    for (int i = 0; i < base_pileup->bases().length(); ++i) {
        if (base_pileup->bases()[i] == '+') {
            ++insert_count;
        }
    }
    if (base_pileup->num_bases() > 0 &&
        (double)insert_count / (double)base_pileup->num_bases() > _max_insert_frac_read_end) {
        base_pileup->set_num_bases(0);
        base_pileup->mutable_bases()->clear();
        base_pileup->mutable_qualities()->clear();
    }
}

Pileups& Pileups::merge(Pileups& other) {
    for (auto& p : other._node_pileups) {
        insert_node_pileup(p.second);
    }
    other._node_pileups.clear();
    for (auto& p : other._edge_pileups) {
        insert_edge_pileup(p.second);
    }
    other._edge_pileups.clear();
    return *this;
}

BasePileup& Pileups::merge_base_pileups(BasePileup& p1, BasePileup& p2) {
    assert(p1.num_bases() == 0 || p2.num_bases() == 0 ||
           p1.ref_base() == p2.ref_base());
    if (p1.num_bases() == 0) {
        p1.set_ref_base(p2.ref_base());
    }
    p1.mutable_bases()->append(p2.bases());
    p1.mutable_qualities()->append(p2.qualities());
    p1.set_num_bases(p1.num_bases() + p2.num_bases());
    p2.set_num_bases(0);
    p2.clear_bases();
    p2.clear_qualities();
    return p1;
}

NodePileup& Pileups::merge_node_pileups(NodePileup& p1, NodePileup& p2) {
    assert(p1.node_id() == p2.node_id());
    for (int i = 0; i < p2.base_pileup_size(); ++i) {
        BasePileup* bp1 = get_create_base_pileup(p1, i);
        BasePileup* bp2 = get_base_pileup(p2, i);
        merge_base_pileups(*bp1, *bp2);
    }
    p2.clear_base_pileup();
    return p1;
}

EdgePileup& Pileups::merge_edge_pileups(EdgePileup& p1, EdgePileup& p2) {
    assert(p1.edge().from() == p2.edge().from());
    assert(p1.edge().to() == p2.edge().to());
    assert(p1.edge().from_start() == p2.edge().from_start());
    assert(p1.edge().to_end() == p2.edge().to_end());
    
    p1.set_num_reads(p1.num_reads() + p2.num_reads());
    p1.mutable_qualities()->append(p2.qualities());
    p2.set_num_reads(0);
    p2.clear_qualities();
    return p1;
}

void Pileups::parse_base_offsets(const BasePileup& bp,
                                 vector<pair<int, int> >& offsets) {
    offsets.clear();
    
    const string& quals = bp.qualities();
    const string& bases = bp.bases();
    char ref_base = ::toupper(bp.ref_base());
    // we can use i to index the quality for the ith row of pileup, but
    // need base_offset to get position of appropriate token in bases string
    int base_offset = 0;
    for (int i = 0; i < bp.num_bases(); ++i) {
        // insert
        if (bases[base_offset] == '+') {
            offsets.push_back(make_pair(base_offset, i < quals.length() ? i : -1));
            int lf = base_offset + 1;
            int rf = lf;
            while (rf < bases.length() && bases[rf] >= '0' && bases[rf] <= '9') {
                ++rf;
            }
            stringstream ss(bases.substr(lf, rf - lf + 1));
            int indel_len;
            ss >> indel_len;
            // ex: +5aaaaa.  rf = lf = 1. indel_len = 5 -> increment 2+0+5=7
            base_offset += 1 + rf - lf + indel_len;
        // delete
        } else if (bases[base_offset] == '-') {
            offsets.push_back(make_pair(base_offset, i < quals.length() ? i : -1));
            int lf = base_offset + 1;
            // eat up four semicolons
            for (int sc_count = 0; sc_count < 4; ++lf) {
                if (bases[lf] == ';') {
                    ++sc_count;
                }
            }
            // and last number
            for (; bases[lf] >= '0' && bases[lf] <= '9'; ++lf);
            base_offset = lf;
        }
        // match / snp
        else {
            offsets.push_back(make_pair(base_offset, i < quals.length() ? i : -1));
            ++base_offset;
        }
    }
    assert(base_offset == bases.length());
}

// transform case of every character in string
void Pileups::casify(string& seq, bool is_reverse) {
    if (is_reverse) {
        transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
    } else {
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    }
}

// make the sam pileup style token
void Pileups::make_match(string& seq, int64_t from_length, bool is_reverse) {
    if (seq.length() == 0) {
        seq = string(from_length, is_reverse ? ',' : '.');
    } else {
        casify(seq, is_reverse);
    }
}

void Pileups::make_insert(string& seq, bool is_reverse) {
    casify(seq, is_reverse);
    stringstream ss;
    ss << "+" << seq.length() << seq; 
    seq = ss.str();
}

void Pileups::make_delete(string& seq, int node_id, int node_offset, bool is_reverse) {
    int delta = is_reverse ? -seq.length() - 1 : seq.length() + 1;
    int dest_offset = node_offset + delta;
    assert(dest_offset >= 0);
    stringstream ss;
    // format : -length;from_start;dest_id;dest_offset;to_end
    ss << "-" << seq.length() << ";" << is_reverse << ";" << node_id << ";"
       << dest_offset << ";" << is_reverse;

    seq = ss.str();
}

void Pileups::append_delete(string& bases, const string& seq) {
    cerr <<"appending " << seq << " to " << bases << endl;
    // this could be streamlined a bit to avoid reparsing
    int d_start = bases.rfind('-');
    assert(d_start >= 0);
    // get existing delete at end of bases
    int old_len, old_offset, old_id;
    bool old_from_start, old_to_end;
    parse_delete(bases.substr(d_start), old_len, old_from_start, old_id, old_offset, old_to_end);

    // parse new delete
    int new_len, new_offset, new_id;
    bool new_from_start, new_to_end;
    parse_delete(seq, new_len, new_from_start, new_id, new_offset, new_to_end);

    // add together
    stringstream merged_del;
    // should be able to use make_delete or something for this...
    merged_del << "-" << (old_len + new_len) << ";" << old_from_start << ";" << new_id << ";"
               << new_offset << ";" << new_to_end;
    
    bases.erase(d_start);
    bases.append(merged_del.str());
}
        
void Pileups::parse_insert(const string& tok, int& len, string& seq, bool& is_reverse) {
    assert(tok[0] == '+');
    int i = 1;
    for (; tok[i] >= '0' && tok[i] <= '9'; ++i);
    stringstream ss;
    ss << tok.substr(1, i - 1);
    ss >> len;
    seq = tok.substr(i, tok.length() - i);
    is_reverse = ::islower(seq[0]);
}

void Pileups::parse_delete(const string& tok, int& len, bool& from_start, int& to_id, int& to_offset,
                           bool& to_end)
{
    assert(tok[0] == '-');
    vector<string> toks;
    split_delims(tok, ";", toks);
    assert(toks.size() == 5);
    len = -atoi(toks[0].c_str());
    from_start = toks[1] != "0";
    to_id = atoi(toks[2].c_str());
    to_offset = atoi(toks[3].c_str());
    to_end = toks[4] != "0";
}
    
bool Pileups::base_equal(char c1, char c2, bool is_reverse) {
    char t1 = ::toupper(c1);
    char t2 = ::toupper(c2);
    return is_reverse ? t1 == reverse_complement(t2) : t1 == t2;
}

char Pileups::extract_match(const BasePileup& bp, int offset) {
    char v = bp.bases()[offset];
    assert(v != '+' && v != '-');
    if (v == ',' || v == '.') {
        return ::toupper(bp.ref_base());
    } else if (::islower(v)) {
        return reverse_complement(::toupper(v));
    }
    return v;
}

// get arbitrary value from offset on forward strand
string Pileups::extract(const BasePileup& bp, int offset) {
    const string& bases = bp.bases();
    if (bases[offset] != '+' && bases[offset] != '-') {
        return string(1, extract_match(bp, offset));
    }
    else if (bases[offset] == '+') {
        string len_str;
        for (int i = offset + 1; bases[i] >= '0' && bases[i] <= '9'; ++i) {
            len_str += bases[i];
        }
        int len = atoi(len_str.c_str());
        // forward strand, return as is
        if (::isupper(bases[offset + 1 + len_str.length()])) {
            return bases.substr(offset, 1 + len_str.length() + len);
        }
        // reverse strand, flip the dna bases part and return upper case
        else {
            string dna = bases.substr(offset + 1 + len_str.length(), len);
            casify(dna, false);
            return string(1, bases[offset]) + len_str + reverse_complement(dna);
        }
    }
    else {
        assert(bases[offset] == '-');
        // todo : consolidate deletion parsing code better than this
        int sc = 0;
        int i = offset;
        for (; sc < 4; ++i) {
            if (bases[i] == ';') {
                ++sc;
            }
        }
        return bases.substr(offset, i - offset + 1);
    }
}

}
