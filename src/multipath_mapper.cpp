//
//  multipath_mapper.cpp
//  
//
//

#include "multipath_mapper.hpp"

namespace vg {
    
    MultipathAligner::MultipathAligner() {
        
    }
    
    MultipathAligner::multipath_align(const Alignment& alignment,
                                      const vector<MaximalExactMatch>& mems){
        
        
        
        // cluster the MEMs
        MultipathClusterer clusterer(alignment, mems, aligner, xgindex, full_length_bonus, max_strand_dist_probes,
                                     max_expected_dist_approx_error);
        vector<vector<pair<MaximalExactMatch* const, pos_t>>> clusters = clusterer.clusters();
        
        // subgraphs around each cluster
        vector<VG*> cluster_graphs;
        cluster_graphs.reserve(clusters.size());
        
        // we will ensure that nodes are in only one cluster, use this to record which one
        unordered_map<id_t, size_t> node_id_to_cluster;
        
        for (size_t i = 0; i < clusters.size(); i++) {
            
            // gather the parameters for subgraph extraction from the MEM hits
            
            vector<pair<MaximalExactMatch* const, pos_t>>& cluster = clusters[i];
            vector<pos_t> positions;
            vector<size_t> forward_max_dist;
            vector<size_t> backward_max_dist;
            
            positions.reserve(cluster.size());
            forward_max_dist.reserve(cluster.size());
            backward_max_dist.reserve(cluster.size());
            
            for (auto mem_hit : cluster) {
                // get the start position of the MEM
                positions.push_back(mem_hit.second);
                // search far enough away to get any hit detectable without soft clipping
                forward_max_dist.push_back(aligner.longest_detectable_gap(alignment, mem_hit.first->end, full_length_bonus)
                                           + (mem_hit.first->end - mem_hit.first->begin));
                backward_max_dist.push_back(aligner.longest_detectable_gap(alignment, mem_hit.first->begin, full_length_bonus));
            }
            
            // TODO: a progressive expansion of the subgraph if the MEM hit is already contained in
            // a cluster graph somewhere?
            
            // extract the subgraph within the search distance
            
            Graph graph;
            algorithms::extract_containing_graph(xgindex, graph, positions, forward_max_dist,
                                                 backward_max_dist, &node_cache);
            
            // check if this subgraph overlaps with any previous subgraph (indicates a probable clustering failure where
            // one cluster was split into multiple clusters)
            unordered_set<size_t> overlapping_graphs;
            for (size_t j = 0; j < graph.node_size(); j++) {
                id_t node_id = graph.node(j).id();
                if (node_id_to_cluster.count(node_id)) {
                    overlapping_graphs.insert(node_id_to_cluster[node_id]);
                }
                else {
                    node_id_to_cluster[node_id] = i;
                }
            }
            
            if (overlapping_graphs.empty()) {
                // there is no overlap with any other graph, suggesting a new unique hit
                
                // hacky socketing into the VG constructor so I can use the move constructor and hopefully
                // avoid reallocating the Graph
                bool passed_graph = false;
                auto pass_graph = [&](Graph& into) {
                    if (passed_graph) {
                        return false;
                    }
                    else {
                        into = std::move(graph);
                        passed_graph = true;
                        return true;
                    }
                };
                cluster_graphs.push_back(new VG(pass_graph));
            }
            else if (overlapping_graphs.size() == 1) {
                // this subgraph overlaps with one other graph, so we merge it in rather than make a new one
                size_t overlap_idx = *overlapping_graphs.begin();
                cluster_graphs[overlap_idx]->extend(graph);
                
                // relabel the cluster of any nodes that were not in the graph we merged into
                Graph& merged_graph = cluster_graphs[overlap_idx].graph;
                for (size_t j = 0; j < merged_graph.node_size(); j++) {
                    node_id_to_cluster[merged_graph.node(j).id()] = overlap_idx;
                }
            }
            else {
                // this subgraph chains together multiple clusters, so we merge all of them into the cluster
                // with the smallest index
                
                size_t min_idx_cluster = std::numeric_limits<size_t>::max();
                for (size_t j : overlapping_graphs) {
                    if (j < min_idx_cluster) {
                        min_idx_cluster = j;
                    }
                }
                
                // merge in the new graph
                cluster_graphs[min_idx_cluster]->extend(graph);
                
                // merge in all the other graphs it connects to and remove them from the list
                overlapping_graphs.erase(min_idx_cluster);
                for (size_t j : overlapping_graphs) {
                    std::swap(cluster_graphs[j], cluster_graphs.back());
                    cluster_graphs[min_idx_cluster]->extend(cluster_graphs.back());
                    delete cluster_graphs.back();
                    cluster_graphs.pop_back();
                }
                
                // relabel the cluster of any nodes that were not in the graph we merged into
                Graph& merged_graph = cluster_graphs[min_idx_cluster].graph;
                for (size_t j = 0; j < merged_graph.node_size(); j++) {
                    node_id_to_cluster[merged_graph.node(j).id()] = min_idx_cluster;
                }
            }
        }
        
        // which MEMs are in play for which cluster?
        unordered_map<VG* vg, vector<pair<MaximalExactMatch* const, pos_t>>> cluster_graph_mems;
        for (const MaximalExactMatch& mem : mems) {
            for (gcsa::node_type hit : mem.nodes) {
                size_t cluster_idx = node_id_to_cluster[id(hit)];
                cluster_graph_mems[cluster_graphs[cluster_idx]].push_back(make_pair(&mem, make_pos_t(hit)));
            }
        }
        
        // record the sequence coverage of each cluster
        unordered_map<VG*, int64_t> cluster_graph_coverage;
        for (VG* vg : cluster_graphs) {
            cluster_graph_coverage[vg] = read_coverage(cluster_graph_mems[vg]);
        }
        
        // sort the cluster graphs descending by unique sequence coverage
        // TODO: figure out relationship between this and the clustering filter
        std::sort(cluster_graphs.begin(), cluster_graphs.end(), [&](const VG* g1, const VG* g2) {
            return cluster_graph_coverage[g1] > cluster_graph_coverage[g2];
        });
        
        // we will preserve paths of up to the maximum length in the graph
        size_t target_length = aligner.longest_detectable_gap(alignment, full_length_bonus) + alignment.sequence().size();
        
        for (VG* vg : cluster_graphs) {
            
            // if we have a cluster graph with small enough MEM coverage compared to the best one
            // we stop producing alternate alignments
            if (cluster_graph_coverage[vg] < mem_coverage_min_ratio * cluster_graph_coverage[cluster_graphs[0]] ) {
                break;
            }
            
//            // pointer to the graph we are going to align against
//            VG* align_graph = vg;
//            
//            // local graphs for unfolding/dagifying as necessary
//            VG unfolded;
//            VG dagified;
//            
//            bool is_acyclic = vg->is_acyclic();
//            bool unrolled = false;
//            unordered_map<id_t, pair<id_t, bool>> trans;
//            
//            if (!is_acyclic) {
//                
//                if (align_graph->has_inverting_edges()) {
//                    unfolded = align_graph->unfold(target_length, trans);
//                    align_graph = &unfolded;
//                    unrolled = true;
//                }
//                
//            }
//            else {
//                
//                unordered_map<id_t, pair<id_t, bool>> unfold_trans;
//                if (align_graph->has_inverting_edges()) {
//                    unfolded = align_graph->unfold(target_length, unfold_trans);
//                    align_graph = &unfolded;
//                    unrolled = true;
//                }
//                
//                dagified = align_graph->dagify(target_length, trans, target_length);
//                align_graph = &dagified;
//                
//                if (!unfold_trans.empty()) {
//                    for (pair<id_t, pair<id_t, bool>>& trans_record : trans) {
//                        pair<id_t, bool>& unfolded_node = unfold_trans[trans_record.second.first];
//                        trans[trans_record.first] = make_pair(unfolded_node.first,
//                                                              trans_record.second.second != unfolded_node.second);
//                    }
//                }
//            }
//            
//            // ensure that graph is in topological order
//            align_graph->sort();
//            
//            // map to project MEMs into the dagified graph
//            unordered_multimap<id_t, pair<id_t, bool>> backward_trans;
//            if (trans.empty()) {
//                const Graph& g = align_graph->graph;
//                for (int64_t i = 0; i < g.node_size(); i++) {
//                    id_t id = g.node(i).id();
//                    backward_trans.insert(make_pair(id, make_pair(g.node(i).id(), false)));
//                }
//            }
//            else {
//                for (pair<id_t, pair<id_t, bool>> trans_record : trans) {
//                    backward_trans.insert(make_pair(trans_record.second.first,
//                                                    make_pair(trans_record.first, trans_record.second.second)));
//                }
//            }
//            
//            vector<pair<MaximalExactMatch* const, pos_t>> forward_hits;
//            vector<pair<MaximalExactMatch* const, pos_t>> reverse_hits;
//            
//            for (const pair<MaximalExactMatch* const, pos_t>& mem_hit : cluster_graph_mems[vg]) {
//                auto range = backward_trans.equal_range(id(mem_hit.second));
//                bool added_forward = false;
//                bool added_reverse = false;
//                for (auto iter = range.first; iter != range.second; iter++) {
//                    if (is_rev(mem_hit.second) == (*iter).second.second && !added_forward) {
//                        forward_hits.push_back(mem_hit);
//                        added_forward = true;
//                    }
//                    else if (is_rev(mem_hit.second) != (*iter).second.second && !added_reverse) {
//                        reverse_hits.push_back(mem_hit);
//                        added_reverse = true;
//                    }
//                }
//            }
//            
//            // sort the hits in descending size order to help filter out partial MEMs
//            
//            std::sort(forward_hits.begin(), forward_hits.end(),
//                      [](const pair<MaximalExactMatch* const, pos_t>& hit_1,
//                         const pair<MaximalExactMatch* const, pos_t>& hit_2) {
//                          return hit_1.first->length() > hit_2.first->length();
//                      });
//            
//            std::sort(reverse_hits.begin(), reverse_hits.end(),
//                      [](const pair<MaximalExactMatch* const, pos_t>& hit_1,
//                         const pair<MaximalExactMatch* const, pos_t>& hit_2) {
//                          return hit_1.first->length() > hit_2.first->length();
//                      });
//            
//            int64_t forward_coverage = read_coverage(forward_hits);
//            int64_t reverse_coverage = read_coverage(reverse_hits);
//            
//            // align to the forward strand if we've detected more hits to it than you would expect to by chance
//            // to a string of equal length to the read (ignoring contiguity) -- this is a permissive standard
//            // but it can be confounded by base errors/polymorphisms that are close enough together that MEMs
//            // are dropped because of the length filter
//            bool align_forward = read_coverage_z_score(forward_coverage, alignment) > z_score_cutoff;
//            // if the graph has been unrolled, we can be more stringent about aligning the reverse strand
//            // because it is likely also represented on the forward strand
//            bool align_backward = unrolled ? reverse_coverage > forward_coverage
//                                           : read_coverage_z_score(reverse_coverage, alignment) > z_score_cutoff;
            
            unordered_map<id_t, pair<id_t, bool> > node_trans;
            VG align_graph = vg->split_strands(node_trans);
            
            unordered_multimap<id_t, pair<id_t, bool> > rev_trans;
            for (const auto& trans_record : node_trans) {
                rev_trans.insert(make_pair(trans_record.second.first,
                                           make_pair(trans_record.first, trans_record.second.second));
            }
            
            vector<pair<MaximalExactMatch* const, pos_t>>& graph_mems = cluster_graph_mems[vg];
            
            // TODO: pick up alignment routine from here
        }
        
        for (VG* vg : cluster_graphs) {
            delete vg;
        }
    }
    
    int64_t MultipathAligner::read_coverage(const vector<pair<MaximalExactMatch* const, pos_t>>& mem_hits) {
        if (mem_hits.empty()) {
            return 0;
        }
        
        vector<pair<string::const_iterator, string::const_iterator>> mem_read_segments;
        mem_read_segments.reserve(mem_hits.size());
        for (auto& mem_hit : mem_hits) {
            mem_read_segments.push_back(mem_hit.first->begin, mem_hit.first->end);
        }
        std::sort(mem_read_segments.begin(), mem_read_segments.end());
        auto curr_begin = mem_read_segments[0].first;
        auto curr_end = mem_read_segments[0].second;
        
        int64_t total = 0;
        for (size_t i = 1; i < trace.size(); i++) {
            if (mem_read_segments[i].first >= curr_end) {
                total += (curr_end - curr_begin);
                curr_begin = mem_read_segments[i].first;
            }
            curr_end = mem_read_segments[i].second;
        }
        return total + (curr_end - curr_begin);
    }
    
    double MultipathAligner::read_coverage_z_score(int64_t coverage, const Alignment& alignment) {
        /* algebraically equivalent to
         *
         *      Coverage - ReadLen / 4
         *  -------------------------------
         *  sqrt(ReadLen * 1/4 * (1 - 1/4))
         * 
         * from the Normal approximation to a Binomal(ReadLen, 1/4)
         */
        double root_len = sqrt(alignment.sequence().size());
        return 0.5773502691896258 * (4.0 * coverage / root_len - root_len);
    }
    
    MultipathAlignment MultipathAligner::set_up_multipath_graph(const Alignment& alignment, VG* vg,
                                                                const vector<pair<MaximalExactMatch* const, pos_t>>& hits,
                                                                const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans,
                                                                const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                                                MultipathAlignmentGraph& multi_aln_graph) {

        
        // map of node ids in the dagified graph to the indices in the matches that contain them
        unordered_map<int64_t, vector<int64_t>> node_matches;
        
        // walk the matches and filter out redundant sub-MEMs
        for (int64_t i = 0; i < hits.size(); i++) {
            
            pair<MaximalExactMatch* const, pos_t>& hit = hits[i];
            
            // the part of the read we're going to match
            string::const_iterator begin = hit.first->begin;
            string::const_iterator end = hit.first->end;
            int64_t mem_length = end - begin;
            // the start of the hit in the original graph
            pos_t& hit_pos = hit.second;
            
            auto hit_range = injection_trans.equal_range(id(hit_pos));
            for (auto iter = hit_range.first; iter != hit_range.second; iter++) {
                // this graph is unrolled/dagified, so all orientations should match
                if ((*iter).second != is_rev(pos)) {
                    continue;
                }
                
                // an id that corresponds to the original node
                id_t injected_id = (*iter).first;
                
                // check all MEMs that traversed this node to see if this is a redundant sub-MEM
                bool is_partial_mem = false;
                for (int64_t j : node_matches[id(hit_pos)]) {
                    ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
                    
                    int64_t relative_offset = begin - match_node.begin;
                    if (relative_offset < 0 || relative_offset + (end - begin) >= match_node.end - match_node.begin) {
                        // the hit does not fall on the same section of the read as the other match, so
                        // it cannot be contained in it
                        continue;
                    }
                    
                    Path& path = match_node.path;
                    
                    // if this is a partial MEM, we should be able to predict its hit location by traversing the path
                    // of the parent MEM by a distance equal to the relative offset
                    
                    int64_t prefix_length = 0;
                    for (size_t k = 0; k < path.mapping_size(); k++) {
                        if (prefix_length > relative_offset) {
                            break;
                        }
                        const Mapping& mapping = path.mapping(k);
                        // the length through this mapping
                        int64_t prefix_through_length = prefix_length + mapping.from_length();
                        if (prefix_through_length > relative_offset) {
                            // we cross the relative offset on this node, so check if the path is in the predicted
                            // position for a redundant sub-MEM
                            id_t node_id_here = mapping.position().node_id();
                            partial_mem = partial_mem ||
                                          (injected_id == node_id_here
                                           && prefix_length + offset(hit_pos) == relative_offset
                                           && projection_trans[node_id_here].second == is_rev(hit_pos));
                        }
                        prefix_length = prefix_through_length;
                    }
                    if (is_partial_mem) {
                        break;
                    }
                }
                
                // don't walk the match of false partial hits
                if (is_partial_mem) {
                    continue;
                }
                
                // stack for DFS, each record contains tuples of (read begin, node offset, next node index, next node ids)
                vector<tuple<string::const_iterator, size_t, size_t, vector<NodeTraversal>>> stack;
                stack.emplace_back(begin, offset(hit_pos), 0,
                                   vector<NodeTraversal>{NodeTraversal(vg->get_node(injected_id))});
                
                while (!stack.empty()) {
                    auto back& = stack.back();
                    if (get<2>(back) == get<3>(back).size()) {
                        stack.pop_back();
                    }
                    NodeTraversal trav = get<3>(back)[get<2>(back)];
                    
                    const string& node_seq = trav.node->sequence();
                    size_t node_idx = get<1>(back);
                    string::const_iterator read_iter = get<0>(back);
                    
                    // look for a match along the entire node sequence
                    for (; node_idx < node_seq.size() && read_iter != end; node_idx++, read_iter++) {
                        if (node_seq[node_idx] != *read_iter) {
                            break;
                        }
                    }
                    
                    if (read_iter == end) {
                        // finished walking match
                        break;
                    }
                    else if (node_idx == node_seq.size()) {
                        // matched entire node
                        auto new_back = stack.emplace_back(read_iter, 0, 0, vector<NodeTraversal>());
                        vg->nodes_next(trav, get<3>(*new_back));
                    }
                    else {
                        // match did not complete node or finish, this is a miss
                        get<2>(back)++;
                    }
                }
                
                if (!stack.empty()) {
                    // we left the trace in the stack, which means we found a complete match
                    int64_t match_node_idx = multi_aln_graph.match_nodes.size();
                    ExactMatchNode& match_node = *(multi_aln_graph.match_nodes.emplace_back());
                    Path& path = match_node.path;
                    match_node.begin = begin;
                    match_node.end = end;
                    
                    // walk out the match
                    int32_t rank = 1;
                    for (auto search_record : stack) {
                        int64_t offset = get<1>(search_record);
                        Node* node = get<3>(search_record)[get<2>(search_record)].node;
                        int64_t length = node->sequence().size() - offset;
                        
                        Mapping* mapping = path.add_mapping();
                        mapping->set_from_length(length);
                        mapping->set_to_length(length);
                        mapping->set_rank(rank);
                        
                        // note: the graph is dagified and unrolled, so all hits should be on the forward strand
                        Position* position = mapping->mutable_position();
                        position->set_node_id(node->id());
                        position->set_offset(offset);
                        
                        // record that each node occurs in this match so we can filter out sub-MEMs
                        node_matches[node->id()].push_back(match_node_idx);
                        
                        rank++;
                    }
                }
            }
        }
        
        // now we calculate reachability between the walked matches so we know which ones
        // to connect with intervening alignments
        
        auto start_offset = [&](size_t idx) {
            return multi_aln_graph.match_nodes[idx].path.mapping(0).position().offset();
        };
        
        auto end_offset = [&](size_t idx) {
            Path& path = multi_aln_graph.match_nodes[idx].path;
            const Mapping& mapping = path.mapping(path.mapping_size() - 1);
            return mapping.position().offset() + mapping_from_length(mapping);
        };
        
        auto start_node_id = [&](size_t idx) {
            return multi_aln_graph.match_nodes[idx].path.mapping(0).position().node_id();
        };
        
        auto end_node_id = [&](size_t idx) {
            Path& path = multi_aln_graph.match_nodes[idx].path;
            return path.mapping(path.mapping_size() - 1).position().node_id();
        };

        // record the start and end node ids of every exact match
        unordered_map<id_t, vector<size_t>> exact_match_starts;
        unordered_map<id_t, vector<size_t>> exact_match_ends;
        for (size_t i = 0; i < multi_aln_graph.match_nodes.size(); i++) {
            Path& path = multi_aln_graph.match_nodes[i].path;
            exact_match_starts[path.mapping(0).position().node_id()].push_back(i);
            exact_match_ends[path.mapping(path.mapping_size() - 1).position().node_id()].push_back(i);
        }
        
        // sort the MEMs starting and ending on each node in node sequence order
        for (pair<id_t, vector<size_t>>& node_match_starts : exact_match_starts) {
            std::sort(node_match_starts.second.begin(), node_match_starts.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                          return start_offset(idx_1) < start_offset(idx_2);
                      });
        }
        for (pair<id_t, vector<size_t>>& node_match_ends : exact_match_ends) {
            std::sort(node_match_ends.second.begin(), node_match_ends.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                          return end_offset(idx_1) < end_offset(idx_2);
                      });
        }
        
        // some structures we will fill out with DP:
        
        // for each node, the starts and ends of MEMs that can reach this node without passing any other
        // MEM starts or ends
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_ends;
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_starts;
        
        // for the start of each MEM, the starts and ends of other MEMs that can reach it without passing any
        // other start or end
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_start;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_start;
        
        // for the end of each MEM, the ends of other MEMs that can reach it without passing any
        // other start or end
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_end;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_end;
        
        // note: graph has been sorted into topological order
        
        Graph& graph = vg->graph;
        for (int64_t i = 0; i < graph.node_size(); i++) {
            Node* node = graph.mutable_node(i);
            id_t node_id = node->id();
            
            size_t node_length = node->sequence().size();
            
            // do any MEMs start or end on this node?
            bool contains_starts = exact_match_starts.count(node_id);
            bool contains_ends = exact_match_ends.count(node_id);
            
            // we will use DP to carry reachability information forward onto the next nodes
            vector<NodeTraversal> nexts;
            vg->nodes_from(NodeTraversal(node), nexts);
            
            if (contains_starts && contains_ends) {
                // since there are both starts and ends on this node, we have to traverse both lists simultaneously
                // to assess reachability within the same node
                
                vector<size_t>& ends = exact_match_starts[node_id];
                vector<size_t>& starts = exact_match_starts[node_id];
                
                // find the range of starts and ends in the list with the same offset
                
                size_t start_range_begin = 0;
                size_t start_range_end = 0;
                size_t end_range_begin = 0;
                size_t end_range_end = 0;
                
                size_t curr_start_offset = start_offset(starts[start_range_begin]);
                size_t curr_end_offset = end_offset(ends[end_range_begin]);
                size_t prev_offset = 0;
                
                bool prev_is_end = curr_end_offset <= curr_start_offset;
                while (end_offset(starts[start_range_end]) == curr_end_offset) {
                    end_range_end++;
                    if (end_range_end == ends.size()) {
                        break;
                    }
                }
                while (start_offset(starts[start_range_end]) == curr_start_offset) {
                    start_range_end++;
                    if (start_range_end == starts.size()) {
                        break;
                    }
                }
                
                // connect the first range of starts or ends to the incoming starts and ends
                
                size_t prev_end_range_begin = end_range_begin;
                size_t prev_start_range_begin = start_range_begin;
                if (prev_is_end) {
                    for (size_t j = end_range_begin; j < end_range_end; j++) {
                        for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                            reachable_ends_from_end[ends[j]].push_back(incoming_end);
                        }
                        for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                            reachable_starts_from_end[ends[j]].push_back(incoming_start);
                        }
                    }
                    
                    end_range_begin = end_range_end;
                    prev_offset = curr_end_offset;
                    curr_end_offset = end_offset(ends[end_range_begin]);
                    if (end_range_begin != ends.size()) {
                        while (end_offset(starts[start_range_end]) == curr_end_offset) {
                            end_range_end++;
                            if (end_range_end == ends.size()) {
                                break;
                            }
                        }
                    }
                }
                else {
                    for (size_t j = start_range_begin; j < start_range_end; j++) {
                        for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                            reachable_ends_from_start[starts[j]].push_back(incoming_end);
                        }
                        for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                            reachable_starts_from_start[starts[j]].push_back(incoming_start);
                        }
                    }
                    
                    start_range_begin = start_range_end;
                    prev_offset = curr_start_offset;
                    curr_start_offset = start_offset(starts[start_range_begin]);
                    if (start_range_begin != starts.size()) {
                        while (start_offset(starts[start_range_end]) == curr_start_offset) {
                            start_range_end++;
                            if (start_range_end == starts.size()) {
                                break;
                            }
                        }
                    }
                }
                
                // iterate along ranges of starts or ends in order of their offsets
                
                while (start_range_begin < starts.size() && end_range_begin < ends.size()) {
                    if (curr_end_offset <= curr_start_offset) {
                        
                        size_t dist_between = curr_end_offset - prev_offset;
                        
                        // connect this range to the previous range
                        if (prev_is_end) {
                            for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                                for (size_t k = end_range_begin; k < end_range_end; k++) {
                                    reachable_ends_from_end[ends[k]].push_back(make_pair(ends[j], dist_between));
                                }
                            }
                        }
                        else {
                            for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                                for (size_t k = end_range_begin; k < end_range_end; k++) {
                                    reachable_starts_from_end[ends[k]].push_back(make_pair(starts[j], dist_between));
                                }
                            }
                        }
                        
                        // record the properties of this range
                        prev_end_range_begin = end_range_begin;
                        prev_is_end = true;
                        
                        // advance to the next range
                        end_range_begin = end_range_end;
                        prev_offset = curr_end_offset;
                        curr_end_offset = end_offset(ends[end_range_begin]);
                        if (end_range_begin != ends.size()) {
                            while (end_offset(starts[start_range_end]) == curr_end_offset) {
                                end_range_end++;
                                if (end_range_end == ends.size()) {
                                    break;
                                }
                            }
                        }
                    }
                    else {
                        
                        size_t dist_between = curr_end_offset - prev_offset;
                        
                        // connect this range to the previous range
                        if (prev_is_end) {
                            for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                                for (size_t k = start_range_begin; k < start_range_end; k++) {
                                    reachable_ends_from_start[start[k]].push_back(make_pair(ends[j], dist_between));
                                }
                            }
                        }
                        else {
                            for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                                for (size_t k = start_range_begin; k < start_range_end; k++) {
                                    reachable_starts_from_start[starts[k]].push_back(make_pair(starts[j], dist_between));
                                }
                            }
                        }
                        
                        // record the properties of this range
                        prev_start_range_begin = start_range_begin;
                        prev_is_end = false;
                        
                        // advance to the next range
                        start_range_begin = start_range_end;
                        prev_offset = curr_start_offset;
                        curr_start_offset = start_offset(starts[start_range_begin]);
                        if (start_range_begin != starts.size()) {
                            while (start_offset(starts[start_range_end]) == curr_start_offset) {
                                start_range_end++;
                                if (start_range_end == starts.size()) {
                                    break;
                                }
                            }
                        }
                    }
                }
                
                // finish off the list of starts on this node
                
                while (start_range_begin < starts.size()) {
                    
                    size_t dist_between = curr_start_offset - prev_offset;
                    
                    if (prev_is_end) {
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = start_range_begin; k < start_range_end; k++) {
                                reachable_ends_from_start[start[k]].push_back(make_pair(ends[j], dist_between));
                            }
                        }
                    }
                    else {
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = start_range_begin; k < start_range_end; k++) {
                                reachable_starts_from_start[starts[k]].push_back(make_pair(starts[j], dist_between));
                            }
                        }
                    }
                    
                    prev_start_range_begin = start_range_begin;
                    start_range_begin = start_range_end;
                    prev_offset = curr_start_offset;
                    curr_start_offset = start_offset(starts[start_range_begin]);
                    prev_is_end = false;
                    
                    if (start_range_begin != starts.size()) {
                        while (start_offset(starts[start_range_end]) == curr_start_offset) {
                            start_range_end++;
                            if (start_range_end == starts.size()) {
                                break;
                            }
                        }
                    }
                }
                
                // finish off the list of ends on this node
                
                while (end_range_begin < ends.size()) {
                    
                    size_t dist_between = curr_start_offset - prev_offset;
                    
                    if (prev_is_end) {
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = end_range_begin; k < end_range_end; k++) {
                                reachable_ends_from_end[ends[k]].push_back(make_pair(ends[j], dist_between));
                            }
                        }
                    }
                    else {
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = end_range_begin; k < end_range_end; k++) {
                                reachable_starts_from_end[ends[k]].push_back(make_pair(starts[j], dist_between));
                            }
                        }
                    }
                    
                    prev_end_range_begin = end_range_begin;
                    end_range_begin = end_range_end;
                    prev_offset = curr_end_offset;
                    curr_end_offset = end_offset(ends[end_range_begin]);
                    prev_is_end = true;
                    
                    if (end_range_begin != ends.size()) {
                        while (end_offset(starts[start_range_end]) == curr_end_offset) {
                            end_range_end++;
                            if (end_range_end == ends.size()) {
                                break;
                            }
                        }
                    }
                }
                
                // carry forward the reachability of the last range onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
                if (prev_is_end) {
                    for (NodeTraversal next : nexts) {
                        unordered_map<size_t, size_t> reachable_ends_next = reachable_ends[next.node->id()];
                        for (size_t j = prev_range_begin; j < ends.size(); j++) {
                            if (reachable_ends_next.count(ends[j])) {
                                reachable_ends_next[ends[j]] = std::min(reachable_ends_next[ends[j]], dist_thru);
                            }
                            else {
                                reachable_ends_next[ends[j]] = dist_thru;
                            }
                        }
                    }
                }
                else {
                    for (NodeTraversal next : nexts) {
                        unordered_map<size_t, size_t> reachable_starts_next = reachable_starts[next.node->id()];
                        for (size_t j = prev_range_begin; j < starts.size(); j++) {
                            if (reachable_starts_next.count(ends[j])) {
                                reachable_starts_next[starts[j]] = std::min(reachable_starts_next[starts[j]], dist_thru);
                            }
                            else {
                                reachable_starts_next[starts[j]] = dist_thru;
                            }
                        }
                    }
                }
            }
            else if (contains_starts) {
                
                // record the incoming starts/ends for the starts on this node
                
                vector<size_t>& starts = exact_match_starts[node_id];
                
                // the starts/ends coming into this node from outside
                size_t range_begin = 0;
                size_t range_end = 0;
                size_t curr_offset = start_offset(starts[range_begin]);
                size_t prev_offset = curr_offset;
                // find the range of starts that are at the first offset
                while (start_offset(starts[range_end]) == curr_offset) {
                    range_end++;
                    if (range_end == starts.size()) {
                        break;
                    }
                }
                // connect the range to the incoming starts/ends
                for (size_t j = range_begin; j < range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                        reachable_starts_from_start[starts[j]].push_back(incoming_start);
                    }
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                        reachable_ends_from_start[starts[j]].push_back(incoming_end);
                    }
                }
                
                // the reachable starts internal to this node
                size_t prev_range_begin = 0;
                range_begin = range_end;
                while (range_begin < starts.size()) {
                    // find the range of starts at this offset
                    prev_offset = curr_offset;
                    curr_offset = start_offset(starts[range_begin]);
                    while (start_offset(starts[range_end]) == curr_offset) {
                        range_end++;
                        if (range_end == starts.size()) {
                            break;
                        }
                    }
                    
                    size_t dist_between = curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    for (size_t j = range_begin; j < range_end; j++) {
                        for (size_t k = prev_range_begin; k < range_begin; k++) {
                            reachable_starts_from_start[starts[j]].push_back(make_pair(starts[k], dist_between));
                        }
                    }
                    prev_range_begin = range_begin;
                    range_begin = range_end;
                }
                
                // this node contains at least one start of a MEM, so carry forward the reachability of all
                // starts at the final offset onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t> reachable_starts_next = reachable_starts[next.node->id()];
                    for (size_t j = prev_range_begin; j < starts.size(); j++) {
                        if (reachable_starts_next.count(ends[j])) {
                            reachable_starts_next[starts[j]] = std::min(reachable_starts_next[starts[j]], dist_thru);
                        }
                        else {
                            reachable_starts_next[starts[j]] = dist_thru;
                        }
                    }
                }
            }
            else if (contains_ends) {
                
                // record the incoming starts/ends for the ends on this node
                
                vector<size_t>& ends = exact_match_ends[node_id];
                
                // the ends coming into this node from outside
                size_t range_begin = 0;
                size_t range_end = 0;
                size_t curr_offset = end_offset(ends[range_begin]);
                size_t prev_offset = curr_offset;
                // find the range of ends that are at the first offset
                while (end_offset(ends[range_end]) == curr_offset) {
                    range_end++;
                    if (range_end == ends.size()) {
                        break;
                    }
                }
                // connect the range to the incoming ends
                for (size_t j = range_begin; j < range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
                        reachable_ends_from_end[ends[j]].push_back(incoming_end);
                    }
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
                        reachable_starts_from_end[ends[j]].push_back(incoming_start);
                    }
                }
                
                // the reachable ends internal to this node
                size_t prev_range_begin = range_begin;
                range_begin = range_end;
                while (range_begin < ends.size()) {
                    // find the range of ends at this offset
                    prev_offset = curr_offset;
                    curr_offset = end_offset(ends[range_begin]);
                    while (end_offset(ends[range_end]) == curr_offset) {
                        range_end++;
                        if (range_end == ends.size()) {
                            break;
                        }
                    }
                    
                    size_t dist_between = curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    for (size_t j = range_begin; j < range_end; j++) {
                        for (size_t k = prev_range_begin; k < range_begin; k++) {
                            reachable_ends_from_end[ends[j]].push_back(make_pair(ends[k], dist_between));
                        }
                    }
                    prev_range_begin = range_begin;
                    range_begin = range_end;
                }
                
                // this node contains at least one end of a MEM, so carry forward the reachability of all
                // ends at the final offset onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t> reachable_ends_next = reachable_ends[next.node->id()];
                    for (size_t j = prev_range_begin; j < ends.size(); j++) {
                        if (reachable_ends_next.count(ends[j])) {
                            reachable_ends_next[ends[j]] = std::min(reachable_ends_next[ends[j]], dist_thru);
                        }
                        else {
                            reachable_ends_next[ends[j]] = dist_thru;
                        }
                    }
                }
            }
            else {
                // this doesn't contain the start or end of any MEM, so we carry forward the reachability
                // into this node onto the next nodes
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t>& reachable_ends_next = reachable_ends[next.node->id()];
                    for (const pair<size_t, size_t>& reachable_end : reachable_ends[node_id]) {
                        size_t dist_thru = reachable_end.second + node_length;
                        if (reachable_ends_next.count(reachable_end.first)) {
                            reachable_ends_next[reachable_end.first] = std::min(reachable_ends_next[reachable_end.first],
                                                                                dist_thru);
                        }
                        else {
                            reachable_ends_next[reachable_end.first] = dist_thru;
                        }
                    }
                    
                    unordered_map<size_t, size_t>& reachable_starts_next = reachable_starts[next.node->id()];
                    for (const pair<size_t, size_t>& reachable_start : reachable_starts[node_id]) {
                        size_t dist_thru = reachable_start.second + node_length;
                        if (reachable_starts_next.count(reachable_start.first)) {
                            reachable_starts_next[reachable_start.first] = std::min(reachable_starts_next[reachable_start.first],
                                                                                    dist_thru);
                        }
                        else {
                            reachable_starts_next[reachable_start.first] = dist_thru;
                        }
                    }
                }
            }
        }
        
        // now we have the reachability information for the start and end of every MEM in the graph. we
        // will use this to navigate between the MEMs in a way that respects graph reachability so that this
        // phase of the algorithm only needs to pay attention to read colinearity and transitive reducibility
        
        vector<unordered_map<size_t, size_t>> noncolinear_shells(multi_aln_graph.match_nodes.size());
        
        unordered_map<size_t, SuffixTree> suffix_trees;
        
        // tuples of (overlap size, index from, index onto)
        vector<tuple<size_t, size_t, size_t>> confirmed_overlaps;
        
        for (size_t i = 0; i < graph.node_size(); i++) {
            id_t node_id = graph.node(i).id();
            
            if (!exact_match_starts.count(node_id)) {
                continue;
            }
            
            for (size_t start : exact_match_starts[node_id]) {
                // traverse all of the reachable starts to find the adjacent ends that might be colinear
                ExactMatchNode& start_node = multi_aln_graph.match_nodes[start];
                
                unordered_map<size_t, size_t>& noncolinear_shell = noncolinear_shells[start];
                
                // pairs of (dist, index)
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> start_queue;
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> end_queue;
                start_queue.emplace(start, 0);
                
                unordered_map<size_t, size_t> traversed_start;
                
                while (!start_queue.empty()) {
                    pair<size_t, size_t> start_here = start_queue.top();
                    start_queue.pop();
                    if (traversed_start.count(start_here.first)) {
                        continue;
                    }
                    traversed_start[start_here.first] = start_here.second;
                    
                    for (const pair<size_t, size_t>& end : reachable_ends_from_start[start_here.second]) {
                        end_queue.emplace(start_here.first + end.second, end.first);
                    }
                    
                    for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here.second]) {
                        start_queue.emplace(start_here.first + start_next.second, start_here.second);
                    }
                }
                
                unordered_map<size_t, size_t> traversed_end;
                
                while (!end_queue.empty()) {
                    size_t candidate_end = end_queue.top().second;
                    size_t candidate_dist = end_queue.top().first;
                    candidate_queue.pop();
                    if (traversed_end.count(candidate_end)) {
                        continue;
                    }
                    traversed_end[candidate_end] = candidate_dist;
                    
                    ExactMatchNode& candidate_end_node = multi_aln_graph.match_nodes[candidate_end];
                    
                    if (candidate_end_node.end <= start_node.begin) {
                        // these MEMs are colinear, so connect them
                        candidate_end_node.edges.push_back(make_pair(start, candidate_dist));
                        
                        // skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                        // this connection
                        for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
                            end_queue.emplace(candidate_dist + shell_pred.second, shell_pred.first);
                        }
                    }
                    else {
                        // could removing an overlap restore colinearity on the read?
                        bool overlap_rescued = false;
                        if (start_node.end > candidate_end.end && candidate_end.begin < start_node.begin) {
                            
                            // make a suffix tree or get a pre-made one
                            if (!suffix_trees.count(candidate_end)) {
                                suffix_trees[candidate_end] = SuffixTree(candidate_end_node.begin, candidate_end_node.end);
                            }
                            SuffixTree& suffix_tree = suffix_trees[candidate_end];
                            
                            // compute the overlap
                            size_t overlap = suffix_tree.longest_overlap(start_node.begin, start_node.end);
                            
                            // is the overlap sufficient to restore colinearity?
                            if (candidate_end_node.end <= start_node.begin + overlap) {
                                confirmed_overlaps.emplace_back(overlap, candidate_end, start);
                                overlap_rescued = true;
                            }
                        }
                        
                        if (!overlap_rescued) {
                            // these MEMs are noncolinear, so add this predecessor to the shell
                            noncolinear_shell.insert(candidate_end);
                            
                            // there is no connection to block further connections back, so any of this MEMs
                            // predecessors could still be colinear
                            
                            // find the ends that can reach it directly
                            unordered_set<size_t> predecessors;
                            for (size_t pred_end : reachable_ends_from_end[candidate_end]) {
                                predecessors.insert(pred_end);
                            }
                            
                            // traverse backward through any starts to find more ends that can reach this MEM
                            unordered_set<size_t> enqueued;
                            list<size_t> queue;
                            
                            // initialize the queue with the immediate start neighbors
                            for (size_t pred_start : reachable_starts_from_end[candidate_end]) {
                                enqueued.insert(pred_start);
                                queue.push_back(pred_start);
                            }
                            
                            // traverse backwards through starts stopping at any ends
                            while (!queue.empty()) {
                                size_t start_here = queue.front();
                                queue.pop_front();
                                
                                for (size_t pred_end : reachable_ends_from_start[start_here]) {
                                    predecessors.insert(pred_end);
                                }
                                for (size_t start_next : reachable_starts_from_start[start_here]) {
                                    if (!enqueued.count(start_here)) {
                                        start_queue.push_back(start_next);
                                        enqueued.insert(start_next);
                                    }
                                }
                            }
                            
                            // add any new predecessors we discovered in the traversal to the candidate queue
                            for (size_t pred : predecessors) {
                                if (!adjacent_ends.count(pred)) {
                                    adjacent_ends.insert(pred);
                                    candidate_queue.push_back(pred);
                                }
                            }
                        }
                    }
                }
                
                // record the node IDs that this MEM's path traverses
                unordered_set<id_t> match_path_nodes;
                Path& match_path = multi_aln_graph.match_nodes[start].path;
                for (size_t j = 0; j < match_path.mapping_size(); j++) {
                    match_path_nodes.insert(path.mapping(j).position().node_id());
                }
                
                // each record in the queue is (start/end idx, is an end?)
                list<pair<size_t, bool>> path_queue{make_pair(start, true)};
                unordered_set<pair<size_t, bool>> path_enqueued{path_queue.front()};
                
                unordered_set<size_t> candidate_overlap_preds;
                
                while (!path_queue.empty()) {
                    // TODO: should I also be adding these to the noncolinear shell?
                    
                    size_t idx = path_queue.front().first;
                    bool at_end = path_queue.front().second;
                    path_queue.pop_front();
                    
                    // stop searching backwards at the start of the MEM
                    if (idx == start && !at_end) {
                        break;
                    }
                    
                    if (at_end) {
                        // this end is on the path of the MEM, it might be overlap-colinear
                        candidate_overlap_preds.insert(idx);
                        
                        // add the start and end predecessors that fall on the path of the MEM to the queue
                        for (size_t start_pred : reachable_starts_from_end[idx]) {
                            if (!path_enqueued.count(make_pair(start_pred, false)) &&
                                match_path_nodes.count(start_node_id(start_pred))) {
                                
                                path_queue.push_back(make_pair(start_pred, false));
                                path_enqueued.insert(path_queue.back());
                            }
                        }
                        for (size_t end_pred : reachable_ends_from_end[idx]) {
                            if (!path_enqueued.count(make_pair(end_pred, true)) &&
                                match_path_nodes.count(end_node_id(end_pred))) {
                                
                                path_queue.push_back(make_pair(end_pred, true));
                                path_enqueued.insert(path_queue.back());
                            }
                        }
                    }
                    else {
                        // this start is on the path of the MEM, it cannot be overlap-colinear even if
                        // the end is also on the path
                        if (candidate_overlap_preds.count(idx)) {
                            candidate_overlap_preds.erase(idx);
                        }
                        
                        // add the start and end predecessors that fall on the path of the MEM to the queue
                        for (size_t start_pred : reachable_starts_from_start[idx]) {
                            if (!path_enqueued.count(make_pair(start_pred, false)) &&
                                match_path_nodes.count(start_node_id(start_pred))) {
                                
                                path_queue.push_back(make_pair(start_pred, false));
                                path_enqueued.insert(path_queue.back());
                            }
                        }
                        for (size_t end_pred : reachable_ends_from_start[idx]) {
                            if (!path_enqueued.count(make_pair(end_pred, true)) &&
                                match_path_nodes.count(end_node_id(end_pred))) {
                                
                                path_queue.push_back(make_pair(end_pred, true));
                                path_enqueued.insert(path_queue.back());
                            }
                        }
                    }
                }
                
                for (size_t overlap_candidate : candidate_overlap_preds) {
                    ExactMatchNode& overlap_node = multi_aln_graph.match_nodes[overlap_candidate];
                    
                    // make a suffix tree or get a pre-made one
                    if (!suffix_trees.count(overlap_candidate)) {
                        suffix_trees[overlap_candidate] = SuffixTree(overlap_node.begin, overlap_node.end);
                    }
                    SuffixTree& suffix_tree = suffix_trees[overlap_candidate];
                    
                    // compute the overlap
                    size_t overlap = suffix_tree.longest_overlap(start_node.begin, start_node.end);
                    
                    // traverse backward down the path to find the trim point for the overlap
                    Path& overlap_path = overlap_node.path;
                    int64_t remaining = overlap;
                    size_t mapping_idx = overlap_path.mapping_size() - 1;
                    int64_t mapping_len = mapping_from_length(overlap_path.mapping(mapping_idx));
                    while (remaining > mapping_len) {
                        remaining -= mapping_len;
                        mapping_idx--;
                        mapping_len = mapping_from_length(overlap_path.mapping(mapping_idx));;
                    }
                    
                    // get the final position of the path prefix
                    id_t final_id = from_path.mapping(mapping_idx).position().node_id();
                    const Position& first_position = onto_node.path.mapping(0).position();
                    
                    //
                    if (!match_path_nodes.count(final_id) ||
                        (first_position.node_id() == final_id && mapping_len - remaining <= first_position.offset())) {
                        
                        confirmed_overlaps.emplace_back(overlap, overlap_candidate, start);
                        
                    }
                }
            }
        }
        
        // sort in descending order of overlap length
        std::sort(confirmed_overlaps.begin(), confirmed_overlaps.end(),
                  std::greater<tuple<size_t, size_t, size_t>>());
        
        
        // split up each node with an overlap edge onto it
        for (auto overlap_record : confirmed_overlaps) {
            ExactMatchNode& from_node = multi_aln_graph.match_nodes[get<1>(overlap_record)];
            ExactMatchNode& onto_node = multi_aln_graph.match_nodes[get<2>(overlap_record)];
            
            size_t suffix_idx = multi_aln_graph.match_nodes.size();
            
            ExactMatchNode& suffix_node = *(multi_aln_graph.match_nodes.emplace_back());
            
            // transfer the outgoing edges onto the new node
            suffix_node.edges = std::move(onto_node.edges);
            
            // clear the old edges and add a single edge to the suffix
            onto_node.edges.clear();
            onto_node.edges.push_back(suffix_idx);
            
            // add the overlap edge
            from_node.edges.push_back(suffix_idx);
            
            // store the full path and remove it from the node
            Path full_path = std::move(onto_node.path);
            onto_node.path.Clear();
            
            // add mappings from the path until reaching the overlap point
            int64_t remaining = get<0>(overlap_record);
            int64_t mapping_idx = 0;
            int64_t mapping_len = mapping_from_length(full_path.mapping(mapping_idx));
            while (remaining >= mapping_len) {
                *onto_node.path.add_mapping() = full_path.mapping(mapping_idx);
                mapping_idx++;
                mapping_len = mapping_from_length(full_path.mapping(mapping_idx));
            }
            
            if (remaining) {
                // the overlap point is in the middle of a node, need to split a mapping
                
                const Mapping& split_mapping = full_path.mapping(mapping_idx);
                
                // add the prefix of the mapping to the original node
                Mapping* prefix_split = onto_node.path.add_mapping();
                prefix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                prefix_split->set_from_length(remaining);
                prefix_split->set_to_length(remaining);
                
                // add the suffix of the mapping to the new node
                Mapping* suffix_split = suffix_node.path.add_mapping();
                suffix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                suffix_split->mutable_position()->set_offset(remaining);
                suffix_split->set_from_length(mapping_len - remaining);
                suffix_split->set_to_length(mapping_len - remaining);
                
                mapping_idx++;
            }
            
            // add the remaining mappings to the suffix node
            for (; mapping_idx < full_path.mapping_size(); mapping_idx) {
                *suffix_node.path.add_mapping() = full_path.mapping(mapping_idx);
            }
        }
        
    }
    
    MultipathClusterer::MultipathClusterer(const Alignment& alignment,
                                           const vector<MaximalExactMatch>& mems,
                                           const BaseAligner& aligner,
                                           xg::XG& xgindex,
                                           int8_t full_length_bonus,
                                           size_t max_strand_dist_probes,
                                           size_t max_expected_dist_approx_error) : aligner(aligner) {
        
        // the maximum graph distances to the right and left sides detectable from each node
        vector<pair<size_t, size_t>> maximum_detectable_gaps;
        
        // there generally will be at least as many nodes as MEMs, so we can speed up the reallocation
        nodes.reserve(mems.size());
        maximum_detectable_gaps.reserve(mems.size());
        
        for (const MaximalExactMatch& mem : mems) {
            
            // calculate the longest gaps we could detect to the left and right of this MEM
            pair<size_t, size_t> max_gaps(aligner.longest_detectable_gap(alignment, mem.begin, full_length_bonus),
                                          aligner.longest_detectable_gap(alignment, mem.end, full_length_bonus));
            
            for (gcsa::node_type mem_hit : mem.nodes) {
                nodes.emplace_back(mem, make_pos_t(mem_hit), mem_score);
                maximum_detectable_gaps.push_back(max_gaps);
            }
        }
        
        // now we use the distance approximation to cluster the MEM hits according to the strand
        // they fall on using the oriented distance estimation function
        
        // initialize with every hit in its own strand cluster
        vector<vector<size_t>> strand_clusters(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++) {
            strand_clusters[i].push_back(i);
        }
        
        // for recording the distance of any pair that we check with a finite distance
        unordered_map<pair<size_t, size_t>, int64_t> recorded_finite_dists;
        // for recording the number of times elements of a strand cluster have been compared
        // and found an infinite distance
        map<pair<size_t, size_t>, size_t> num_infinite_dists;
        
        // when we merge strand clusters, we add a pointer to the strand they merged into
        // essentially, this become trees with pointers upward that we can traverse to
        // find the unmerged strand cluster any hit currently belongs to
        unordered_map<size_t, size_t> merged_strand;
        auto get_merged_strand = [&merged_strand](size_t i) {
            vector<size_t> path;
            while (merged_strand.count(i)) {
                path.push_back(i);
                i = merged_strand[i];
            }
            // path compression (similar to union find)
            for (size_t j : path) {
                merged_strand[j] = i;
            }
            return i;
        }
        
        // make a randomly shuffled list of pairs of nodes to compare distances in
        // TODO: is there a way to do this without generating all pairs?
        vector<pair<size_t, size_t>> comparisons;
        comparisons.reserve((nodes.size() * (nodes.size() - 1)) / 2);
        for (int64_t i = 1; i < nodes.size(); i++) {
            for (int64_t j = 0; j < i; j++) {
                comparisons.emplace_back(j, i):
            }
        }
        std::default_random_engine gen(std::random_device());
        std::shuffle(comparisons.begin(), comparisons.end(), gen);
        
        for (pair<size_t, size_t>& node_pair : comparisons) {
            
            // TODO: I don't love this component, it makes the algorithm O(n^2 log n) in number of MEM hits
            size_t strand_1 = get_merged_strand(node_pair.first);
            size_t strand_2 = get_merged_strand(node_pair.second);

            if (strand_1 == strand_2) {
                // these are already identified as on the same strand, don't need to do it again
                continue;
            }
            
            if (num_infinite_dists[make_pair(strand_1, strand_2)] >= max_strand_dist_probes) {
                // we've already checked multiple distances between these strand clusters and
                // none have returned a finite distance, so we conclude that they are in fact
                // on separate clusters and decline to check any more distances
                continue;
            }
            
            pos_t& pos_1 = nodes[node_pair.first];
            pos_t& pos_2 = nodes[node_pair.second];
            
            int64_t oriented_dist = xg_index.closest_shared_path_oriented_distance(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                                                   id(pos_2), offset(pos_2), is_rev(pos_2));
            
            if (oriented_dist == std::numeric_limits<int64_t>::max()) {
                // distance is estimated at infinity, so these are either on different strands
                // or the path heuristic failed to find a shared path
                
                num_infinite_dists[make_pair(strand_1, strand_2)]++;
                num_infinite_dists[make_pair(strand_2, strand_1)]++;
            }
            else {
                // the distance is finite, so merge the strand clusters
                
                recorded_finite_dists[node_pair] = oriented_dist;
                
                // add the smaller cluster onto the larger one to minimize copies
                vector<size_t>* smaller_clust, larger_clust;
                
                // ensure that the first node in the pair belongs to the larger cluster
                if (strand_clusters[node_pair.first].size() < strand_clusters[node_pair.second].size()) {
                    std::swap(node_pair.first, node_pair.second);
                    std::swap(strand_1, strand_2);
                }
                
                auto& clust_1 = strand_clusters[node_pair.first];
                auto& clust_2 = strand_clusters[node_pair.second];
                
                clust_1.insert(clust_1.end(), clust_2.begin(), clust_2.end());
                
                // choose one of the strand clusters at random to remove (makes the merged_strand
                // tree have expected height O(log n))
                size_t strand_retaining, strand_removing;
                if (gen() % 2) {
                    merged_strand[node_pair.second] = node_pair.first;
                    std::swap(clust_2, strand_clusters.back());
                    strand_removing = strand_2;
                    strand_retaining = strand_1;
                }
                else {
                    merged_strand[node_pair.first] = node_pair.second;
                    clust_2 = std::move(clust_1);
                    std::swap(clust_1, strand_clusters.back());
                    strand_removing = strand_1;
                    strand_retaining = strand_2;
                }
                strand_clusters.pop_back();
                
                // collect the number of times this strand cluster has had an infinite distance to other strand
                vector<pair<size_t, size_t>> inf_counts;
                auto end = num_infinite_dists.upper_bound(make_pair(strand_removing, numeric_limits<size_t>::max()));
                auto begin = num_infinite_dists.lower_bound(make_pair(strand_removing, 0));
                for (auto iter = begin; iter != end; iter++) {
                    inf_counts.push_back((*iter).first.second, (*iter).second);
                }
                
                // add these counts to the other strand cluster and remove this one from the map
                for (const pair<size_t, size_t> inf_count : inf_counts) {
                    num_infinite_dists[make_pair(strand_retaining, inf_count.first)] += inf_count.second;
                    num_infinite_dists.erase(make_pair(strand_removing, inf_count.first));
                    num_infinite_dists.erase(make_pair(inf_count.first, strand_removing));
                }
            }
        }
        
        // build the graph of relative distances in adjacency list representation
        // by construction each strand cluster will be an undirected, unrooted tree
        unordered_map<size_t, vector<size_t>> strand_distance_tree;
        for (const auto& dist_record : recorded_finite_dists) {
            strand_distance_tree[dist_record.first.first].push_back(dist_record.first.second);
            strand_distance_tree[dist_record.first.second].push_back(dist_record.first.first);
        }
        
        // now approximate the relative positions along the strand by traversing each tree and
        // treating the distances we estimated as transitive
        vector<unordered_map<size_t, int64_t>> strand_relative_position;
        vector<bool> processed(nodes.size());
        for (const auto& adjacency_record : strand_distance_tree) {
            if (processed[adjacency_record.first]) {
                continue;
            }
            
            size_t strand = get_merged_strand(adjacency_record.first);
            unordered_map<size_t, int64_t>& relative_pos = strand_relative_position[strand];
            
            // arbitrarily make this node the 0 point
            relative_pos[adjacency_record.first] = 0;
            processed[adjacency_record.first] = true;
            
            // traverse the strand's tree with DFS
            list<size_t> queue{adjacency_record.first};
            while (!queue.empty()) {
                size_t curr = queue.back();
                queue.pop_back();
                
                int64_t curr_pos = relative_pos[curr];
                
                for (size_t next : strand_distance_tree[curr]) {
                    if (processed[next]) {
                        continue;
                    }
                    
                    // invert the sign of the distance if we originally measured it in the other order
                    int64_t dist = recorded_finite_dists.count(make_pair(curr, next)) ?
                                   recorded_finite_dists[make_pair(curr, next)] :
                                   -recorded_finite_dists[make_pair(next, furr)];
                    
                    // find the position relative to the previous node we just traversed
                    relative_pos[next] = curr_pos + dist;
                    processed[next] = true;
                    
                    queue.push_back(next);
                }
            }
        }
        
        // now we use the strand clusters and the estimated distances to make the DAG for the
        // approximate MEM alignment
        
        int64_t allowance = max_expected_dist_approx_error;
        for (const unordered_map<size_t, int64_t>& relative_pos : strand_relative_position) {
            
            // sort the nodes by relative position
            vector<pair<int64_t, size_t>> sorted_pos;
            for (const pair<size_t, int64_t>& pos_record : relative_pos) {
                sorted_pos.emplace_back(pos_record.second, pos_record.first);
            }
            std::sort(sorted_pos.begin(), sorted_pos.end());
            
            // find edges within each strand cluster by first identifying the interval of MEMs that meets
            // the graph distance constrant for each MEM and then checking for read colinearity and the
            // reverse distance constraint
            int64_t last_idx = sorted_pos.size() - 1;
            int64_t low = 0, hi = last_idx;
            for (int64_t i = 0; i < sorted_pos.size(); i++) {
                
                int64_t strand_pos = sorted_pos[i].first;
                size_t pivot_idx = sorted_pos[i].second;
                MPCNode& pivot = nodes[pivot_idx];
                int64_t pivot_length = pivot.mem->end - pivot.mem->begin;
                
                // the limits of how far away we might detect edges to add to the clustering graph
                int64_t target_low_pos = strand_pos - allowance;
                int64_t target_hi_pos = strand_pos + pivot_length + maximum_detectable_gaps[pivot_idx].second + allowance;
                
                // move the lower boundary of the search interval to the lowest value inside the
                // the target interval
                while (sorted_pos[low].first < target_low_pos) {
                    low++;
                }
                
                // move the upper boundary of the search interval to the highest value inside the
                // the target interval (this one can move in either direction because the maximum
                // detectable gap changes)
                if (sorted_pos[hi].first > target_hi_pos) {
                    while (sorted_pos[hi].first > target_hi_pos) {
                        hi--;
                    }
                }
                else {
                    while (hi == last_idx ? false : sorted_pos[hi + 1].first <= target_hi_pos) {
                        hi++;
                    }
                }
                
                for (int64_t j = low; j <= hi; j++) {
                    int64_t next_idx = sorted_pos[j].second;
                    MPCNode& next = nodes[next_idx];
                    
                    if (next.mem->begin <= pivot.mem->begin || next.mem->end <= pivot.mem->end) {
                        // the MEMs cannot be colinear along the read (also filters out j == i)
                        continue;
                    }
                    
                    // the length of the sequence in between the MEMs (can be negative if they overlap)
                    int64_t between_length = next.mem->begin - pivot.mem->end;
                    // the estimated distance between the end of the pivot and the start of the next MEM in the graph
                    int64_t graph_dist = max(0, sorted_pos[j].first - strand_pos - pivot_length);
                    // the discrepancy between the graph distance and the read distance
                    int64_t gap_length = abs(graph_dist - between_length);
                    
                    if (gap_length > maximum_detectable_gaps[next_idx].first + allowance) {
                        // the gap between the MEMs is too long to be believable from the next node
                        continue;
                    }
                    
                    int32_t edge_score;
                    if (between_length < 0) {
                        // the MEMs overlap, but this can occur in some insertions and deletions
                        // because the SMEM algorithm is "greedy" in taking up as much of the read
                        // as possible
                        // we can check if this happened with the SuffixTree, but it's expensive
                        // so for now we just give it the benefit of the doubt but adjust the edge
                        // score so that the matches don't get double counted
                        
                        int64_t extra_dist = max(0, gap_length);
                        
                        edge_score = -aligner.match * between_length
                                     + (extra_dist ? -(extra_dist - 1) * aligner.gap_extension - aligner.gap_open : 0);
                    }
                    else if (between_length > graph_dist) {
                        // the read length in between the MEMs is longer than the distance, suggesting a read insert
                        edge_score = -aligner.mismatch * graph_dist - (gap_length - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else if (between_length < graph_dist) {
                        // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                        edge_score = -aligner.mismatch * between_length - (gap_length - 1) * aligner.gap_extension
                                     - aligner.gap_open;
                    }
                    else {
                        // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                        edge_score = -aligner.mismatch * between_length;
                    }
                    
                    // add the edges in
                    pivot.edges_from.emplace_back(next_idx, edge_score);
                    next.edges_to.emplace_back(pivot_idx, edge_score);
                }
            }
        }
    }

    void MultipathClusterer::topological_order(vector<size_t>& order_out) {
        
        // initialize return value
        order_out.clear();
        order_out.resize(nodes.size());
        size_t order_idx = nodes.size() - 1;
        
        // initialize iteration structures
        vector<bool> enqueued = vector<bool>(nodes.size());
        vector<size_t> edge_index = vector<size_t>(nodes.size());
        vector<size_t> stack;
        
        // iterate through starting nodes
        for (size_t init_node_idx = 0; init_node_idx < nodes.size(); init_node_idx++) {
            if (enqueued[init_node_idx]) {
                continue;
            }
            // navigate through graph with DFS
            stack.push_back(init_node_idx);
            enqueued[init_node_idx] = true;
            while (!stack.empty()) {
                size_t node_idx = stack.back();
                if (edge_index[node_idx] < nodes[node_idx].edges_from.size()) {
                    int64_t target_idx = nodes[node_idx].edges_from[node_idx].to_idx;
                    if (enqueued[target_idx]) {
                        edge_index[node_idx]++;
                    }
                    else {
                        stack.push_back(target_idx);
                        enqueued[target_idx] = true;
                    }
                }
                else {
                    // add to topological order in reverse finishing order
                    stack.pop_back();
                    order_out[order_idx] = node_idx;
                    order_idx--;
                }
            }
        }
    }
    
    void MultipathMEMAligner::identify_sources_and_sinks(vector<size_t>& sources_out,
                                                         vector<size_t>& sinks_out) {
        
        sources_out.clear();
        sinks_out.clear();
        
        vector<bool> is_source(nodes.size(), true);
        
        for (size_t i = 0; i < nodes.size(); i++) {
            if (nodes[i].edges_from.empty()) {
                sinks_out.push_back(i);
            }
            
            for (MultipathMEMEdge& edge : nodes[i].edges_from) {
                is_source[edge.to_idx] = false;
            }
        }
        
        for (size_t i = 0; i < nodes.size(); i++) {
            if (is_source[i]) {
                sources_out.push_back(i);
            }
        }
    }
    
    
    
    void MultipathClusterer::connected_components(vector<vector<size_t>>& components_out) {
        
        components_out.clear();
        vector<bool> enqueued(nodes.size());
        
        // check each node in turn to find new components
        for (size_t dfs_start_idx = 0; dfs_start_idx < nodes.size(); dfs_start_idx++) {
            if (enqueued[dfs_start_idx]) {
                // we've already found this node from some component
                continue;
            }
            
            // this node belongs to a component we haven't found yet, use DFS to find the rest
            vector<size_t> stack {dfs_start_idx};
            enqueued[dfs_start_idx] = true;
            components_out.emplace_back(1, dfs_start_idx);
            
            while (!stack.empty()) {
                
                MultipathMEMNode& node = nodes[stack.back()];
                stack.pop_back();
                
                // search in both forward and backward directions
                
                for (MultipathMEMEdge& edge : node.edges_from) {
                    
                    if (!enqueued[edge.to_idx]) {
                        stack.push_back(edge.to_idx);
                        enqueued[edge.to_idx] = true;
                        components_out.back().push_back(edge.to_idx);
                    }
                }
                
                for (MultipathMEMEdge& edge : node.edges_to) {
                    
                    if (!enqueued[edge.to_idx]) {
                        stack.push_back(edge.to_idx);
                        enqueued[edge.to_idx] = true;
                        components_out.back().push_back(edge.to_idx);
                    }
                }
            }
        }
    }
    
    void MultipathClusterer::perform_dp() {
        
        // as in local alignment, minimum score is the score of node itself
        for (size_t i = 0; i < nodes.size(); i++) {
            nodes[i].dp_score = nodes[i].score;
        }
        
        vector<size_t> order;
        topological_order(order);
        
        for (size_t i : order) {
            MPCNode& node = nodes[i];
            
            // for each edge out of this node
            for (MPCEdge& edge : node.edges_from) {
                
                // check if the path through the node out of this edge increase score of target node
                MPCNode& target_node = nodes[edge.to_idx];
                int32_t extend_score = node.score + edge.weight + target_node.score;
                if (extend_score > target_node.dp_score) {
                    target_node.dp_score = extend_score;
                }
            }
        }
    }
    
    vector<vector<pair<MaximalExactMatch* const, pos_t>>> MultipathClusterer::clusters(int32_t max_qual_score) {
        
        vector<vector<MaximalExactMatch* const, pos_t>> to_return;
        if (nodes.size() == 0) {
            // this should only happen if we have filtered out all MEMs, so there are none to cluster
            return to_return;
        }
        
        perform_dp();
        
        // find the weakly connected components, which should correspond to mappings
        vector<vector<size_t>> components;
        connected_components(components);
        
        // find the node with the highest DP score in each connected component
        vector<pair<int32_t, size_t>> component_traceback_ends(components.size(),
                                                               pair<int32_t, size_t>(numeric_limits<int32_t>::min(), 0));
        for (size_t i = 0; i < components.size(); i++) {
            vector<size_t>& component = components[i];
            pair<int32_t, size_t>& traceback_end = component_traceback_ends[i];
            for (size_t j = 0; j < component.size(); j++) {
                int32_t dp_score = nodes[component[i]].dp_score
                if (dp_score > traceback_end.first) {
                    traceback_end.first = dp_score;
                    traceback_end.second = j;
                }
            }
        }
        
        // sort indices in descending order by their highest traceback score
        vector<size_t> order = range_vector(0, components.size());
        std::sort(order.begin(), order.end() [&](const size_t i, const size_t j) {
            return component_traceback_ends[i].first > component_traceback_ends[j].first;
        });
        
        int32_t top_score = component_traceback_ends[order[0]].first;
        
        for (size_t i : order) {
            // get the component and the traceback end
            vector<size_t>& component = components[i];
            size_t trace_idx = component[component_traceback_ends[i].second];
            
            // if this cluster does not look like it even affect the mapping quality of the top scoring
            // cluster, don't bother forming it
            // TODO: this approximation could break down sometimes, need to look into it
            // TODO: is there a way to make the aligner do this? I don't like having this functionality outside of it
            if (4.3429448190325183 * aligner.log_base * (top_score - traceback_end.first) > max_qual_score ) {
                continue;
            }
            
            // traceback until hitting a node that has its own score (indicates beginning of a local alignment)
            vector<size_t> trace{trace_idx};
            while (nodes[trace_idx].dp_score > nodes[trace_idx].score) {
                for (MPCEdge& edge : nodes[trace_idx].edges_to) {
                    if (nodes[edge.to_idx].dp_score + edge.weight == nodes[trace_idx].dp_score) {
                        trace_idx = edge.to_idx;
                        trace.push_back(trace_idx);
                    }
                }
            }
            
            // make a cluster
            to_return.emplace_back();
            auto& cluster = to_return.back();
            for (auto iter = trace.rbegin(); iter != trace.rend(); iter++) {
                MPCNode& node = nodes[*iter];
                cluster.emplace_back(node.mem, node.start_pos);
            }
        }
        
        return to_return;
    }

}


