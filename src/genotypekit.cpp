#include "genotypekit.hpp"

namespace vg {

using namespace std;


CactusSiteFinder::CactusSiteFinder(VG& graph, const string& hint_path_name): graph(graph), hint_path_name(hint_path_name) {
    // Make sure the graph is sorted.
    // cactus needs the nodes to be sorted in order to find a source and sink.
    graph.sort();
}

void CactusSiteFinder::for_each_site_parallel(const function<void(const NestedSite)>& lambda) {

    // Set up our output vector
    vector<NestedSite> to_return;
    
    // get endpoints using node ranks
    pair<NodeSide, NodeSide> source_sink = graph.paths.has_path(hint_path_name) ? 
        get_cactus_source_sink(graph, hint_path_name)
        : get_cactus_source_sink(graph);
        
    // Don't keep going if we can't find sources/sinks
    assert(graph.has_node(source_sink.first.node));
    assert(graph.has_node(source_sink.second.node));

    // todo: use deomposition instead of converting tree into flat structure
    BubbleTree bubble_tree = cactusbubble_tree(graph, source_sink);

    // copy nodes up to bubbles that contain their bubble
    bubble_up_bubbles(bubble_tree);

    bubble_tree.for_each_preorder([&](BubbleTree::Node* node) {
            Bubble& bubble = node->v;
            // cut root to be consistent with superbubbles()
            if (bubble.start != bubble_tree.root->v.start ||
                bubble.end != bubble_tree.root->v.end) {
                set<id_t> nodes{bubble.contents.begin(), bubble.contents.end()};
                NodeTraversal start(graph.get_node(bubble.start.node), !bubble.start.is_end);
                NodeTraversal end(graph.get_node(bubble.end.node), bubble.end.is_end);
                // Fill in a Site. Make sure to preserve original endpoint
                // ordering, because swapping them without flipping their
                // orientation flags will make an inside-out site.
                NestedSite site;
                site.start = start;
                site.end = end;
                for(id_t id : nodes) {
                    site.nodes.insert(graph.get_node(id));
                }
                
                // OpenMP doesn't like to let us use the reference, even though
                // we know it will survive the tasks. We grab a pointer to make
                // it happy.
                auto* lambda_pointer = &lambda;
                
                // Operate on the site
                #pragma omp task
                {
                    (*lambda_pointer)(site);
                }
                
            }
    });
        
    // Don't return until all the lambda tasks are done.
    #pragma omp taskwait    
    
}

double FixedGenotypePriorCalculator::calculate_log_prior(const Genotype& genotype) {
    // Are all the alleles the same?
    bool all_same = true;
    // What is the common allele number (-1 for unset)
    int allele_value = -1;
    for(size_t i = 0; i < genotype.allele_size(); i++) {
        // For each allele in the genotype
        if(allele_value == -1) {
            // On the first one, grab it. Everyone else will have to match.
            allele_value = genotype.allele(i);
        }
        
        if(allele_value != genotype.allele(i)) {
            // There are two distinct allele values in here
            all_same = false;
            break;
        }
    }
    
    // Return the appropriate prior depending on whether the alleles are all the
    // same (homozygous) or not (heterozygous).
    return all_same ? homozygous_prior_ln : heterozygous_prior_ln;
}

}