//
//  snarls.hpp
//
//  Contains object to own Snarls and keep track of their tree relationships as well as utility
//  functions that interact with snarls.
//

#ifndef snarls_hpp
#define snarls_hpp

#include <cstdint>
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <deque>
#include "stream.hpp"
#include "vg.hpp"
#include "handle.hpp"
#include "vg.pb.h"
#include "hash_map.hpp"
#include "cactus.hpp"

using namespace std;

namespace vg {

class SnarlManager;

/**
 * Represents a strategy for finding (nested) sites in a vg graph that can be described
 * by snarls. Polymorphic base class/interface.
 */
class SnarlFinder {
public:
    virtual ~SnarlFinder() = default;
    
    /**
     * Run a function on all root-level NestedSites in parallel. Site trees are
     * passed by value so they have a clear place to live during parallel
     * operations.
     */
    virtual SnarlManager find_snarls() = 0;
};

/**
 * Class for finding all snarls using the base-level Cactus snarl decomposition
 * interface.
 */
class CactusSnarlFinder : public SnarlFinder {
    
    /// Holds the vg graph we are looking for sites in.
    VG& graph;
    
    /// Holds the names of reference path hints
    unordered_set<string> hint_paths;
    
    /// Create a snarl in the given SnarlManager with the given start and end,
    /// containing the given child snarls in the list of chains of children and
    /// the given list of unary children. Recursively creates snarls in the
    /// SnarlManager for the children. Returns a pointer to the finished snarl
    /// in the SnarlManager. Start and end may be empty visits, in which case no
    /// snarl is created, all the child chains are added as root chains, and
    /// null is returned. If parent_start and parent_end are empty Visits, no
    /// parent() is added to the produced snarl.
    const Snarl* recursively_emit_snarls(const Visit& start, const Visit& end,
        const Visit& parent_start, const Visit& parent_end,
        stList* chains_list, stList* unary_snarls_list, SnarlManager& destination);
    
public:
    /**
     * Make a new CactusSnarlFinder to find snarls in the given graph.
     * We can't filter trivial bubbles because that would break our chains.
     *
     * Optionally takes a hint path name.
     */
    CactusSnarlFinder(VG& graph);
    
    /**
     * Make a new CactusSnarlFinder with a single hinted path to base the
     * decomposition on.
     */
    CactusSnarlFinder(VG& graph, const string& hint_path);
    
    /**
     * Find all the snarls with Cactus, and put them into a SnarlManager.
     */
    virtual SnarlManager find_snarls();
    
};

/**
 * Snarls are defined at the Protobuf level, but here is how we define
 * chains as real objects.
 *
 * The SnarlManager is going to have one official copy of each chain stored,
 * and it will give you a pointer to it on demand.
 */
using Chain = vector<const Snarl*>;
    
/**
 * Return true if the first snarl in the given chain is backward relative to the chain.
 */
bool start_backward(const Chain& chain);
    
/**
 * Return true if the last snarl in the given chain is backward relative to the chain.
 */
bool end_backward(const Chain& chain);
    
/**
 * Get the inward-facing start Visit for a chain.
 */
Visit get_start_of(const Chain& chain);
    
/**
 * Get the outward-facing end Visit for a chain.
 */
Visit get_end_of(const Chain& chain);
    
/**
 * We want to be able to loop over a chain and get iterators to pairs of the
 * snarl and its orientation in the chain. So we define some iterators.
 */
struct ChainIterator {
    /// Advance the iterator
    ChainIterator& operator++();
    /// Get the snarl we're at and whether it is backward 
    pair<const Snarl*, bool> operator*() const;
    /// Get a pointer to the thing we get when we dereference the iterator
    const pair<const Snarl*, bool>* operator->() const;
        
    /// We need to define comparison because C++ doesn't give it to us for free.
    bool operator==(const ChainIterator& other) const;
    bool operator!=(const ChainIterator& other) const;
        
    /// Are we a reverse iterator or not?
    bool go_left;
    /// Is the snarl we are at backward or forward in its chain?
    bool backward;
        
    /// What position in the underlying vector are we in?
    Chain::const_iterator pos;
        
    /// What are the bounds of that underlying vector?
    Chain::const_iterator chain_start;
    Chain::const_iterator chain_end;
        
    /// Since we're using backing random access itarators to provide reverse
    /// iterators, we need a flag to see if we are rend (i.e. before the
    /// beginning)
    bool is_rend;
        
    /// When dereferencing, should we flip snarl orientations form the
    /// orientations they appear at in the chain when read left to right?
    bool complement;
        
    /// In order to dereference to a pair with -> we need a place to put the pair so we can have a pointer to it.
    /// Gets lazily set to wherever the iterator is pointing when we do ->
    mutable pair<const Snarl*, bool> scratch;
};
    
/**
 * We define free functions for getting iterators forward and backward through chains.
 */
ChainIterator chain_begin(const Chain& chain);
ChainIterator chain_end(const Chain& chain);
ChainIterator chain_rbegin(const Chain& chain);
ChainIterator chain_rend(const Chain& chain);
    
/// We also define some reverse complement iterators, which go from right to
/// left through the chains, but give us the reverse view. For ecample, if
/// all the snarls are oriented forward in the chain, we will iterate
/// through the snarls in reverse order, with each individual snarl also
/// reversed.
ChainIterator chain_rcbegin(const Chain& chain);
ChainIterator chain_rcend(const Chain& chain);
    
/// We also define a function for getting the ChainIterator (forward or
/// reverse complement) for a chain starting with a given snarl in the given
/// inward orientation
ChainIterator chain_begin_from(const Chain& chain, const Snarl* start_snarl, bool snarl_orientation);
/// And the end iterator for the chain (forward or reverse complement)
/// viewed from a given snarl in the given inward orientation
ChainIterator chain_end_from(const Chain& chain, const Snarl* start_snarl, bool snarl_orientation);
    
/**
 * Allow traversing a graph of nodes and child snarl chains within a snarl
 * within another HandleGraph. Uses its own internal child index because
 * it's used in the construction of snarls to feed to SnarlManagers.
 *
 * Assumes that the chains we get from Cactus are in a consistent order, so
 * the start of the first snarl is the very first thing in the chain, and
 * the end of the last snarl is the very last.
 *
 * We adapt the handle graph abstraction as follows:
 *
 * A chain becomes a single node with the ID and local forward orientation
 * of its first snarl's start.
 *
 * A chain node connects on its left to everything connected to its first
 * start and on its right to everything connected to its last end.
 *
 * A unary snarl becomes a single node, too. It is identified by its
 * boundary node's ID.
 *
 * If you're not using internal connectivity, a chain node or a unary snarl
 * node behaves just like an ordinary node.
 *
 * If you are using internal connectivity, edges are slightly faked:
 *
 * A chain node also sees out its right everything that is out its left if
 * it has a left-left connected snarl before any disconnected snarl.
 *
 * And similarly for the mirror case.
 *
 * All the edges on either side of a unary snarl node are the same.
 *
 * In this part of the code we talk about "heads" (the inward-facing base
 * graph handles used to represent child snarls/chains), and "tails" (the
 * inward-facing ending handles of child chains).
 * 
 */
class NetGraph : public HandleGraph {
public:
        
    /// Make a new NetGraph for the given snarl in the given backing graph,
    /// using the given chains as child chains. Unary snarls are stored as
    /// single-snarl chains just like other trivial chains.
    template<typename ChainContainer>
    NetGraph(const Visit& start, const Visit& end,
             const ChainContainer& child_chains_mixed,
             const HandleGraph* graph,
             bool use_internal_connectivity = false) : NetGraph(start, end, graph, use_internal_connectivity) {
            
        // All we need to do is index the children. They come mixed as real chains and unary snarls.
            
        for (auto& chain : child_chains_mixed) {
            if (chain.size() == 1 && chain.front()->type() == UNARY) {
                // This is a unary snarl wrapped in a chain
                add_unary_child(chain.front());
            } else {
                // This is a real (but possibly trivial) chain
                add_chain_child(chain);
            }
        }
            
    }
        
    /// Make a net graph from the given chains and unary snarls (as pointers) in the given backing graph.
    template<typename ChainContainer, typename SnarlContainer>
    NetGraph(const Visit& start, const Visit& end,
             const ChainContainer& child_chains,
             const SnarlContainer& child_unary_snarls, const HandleGraph* graph,
             bool use_internal_connectivity = false) : NetGraph(start, end, graph, use_internal_connectivity) {
            
        // All we need to do is index the children.
            
        for (const Snarl* unary : child_unary_snarls) {
            add_unary_child(unary);
        }
            
        for (auto& chain : child_chains) {
            add_chain_child(chain);
        }
    }
            
    /// Make a net graph from the given chains and unary snarls (as raw values) in the given backing graph.
    /// Mostly for testing.
    NetGraph(const Visit& start, const Visit& end,
             const vector<vector<Snarl>>& child_chains,
             const vector<Snarl>& child_unary_snarls,
             const HandleGraph* graph,
             bool use_internal_connectivity = false);
            
    /// Look up the handle for the node with the given ID in the given orientation
    virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
    // Copy over the visit version which would otherwise be shadowed.
    using HandleGraph::get_handle;
        
    /// Get the ID from a handle
    virtual id_t get_id(const handle_t& handle) const;
        
    /// Get the orientation of a handle
    virtual bool get_is_reverse(const handle_t& handle) const;
        
    /// Invert the orientation of a handle (potentially without getting its ID)
    virtual handle_t flip(const handle_t& handle) const;
        
    /// Get the length of a node
    virtual size_t get_length(const handle_t& handle) const;
        
    /// Get the sequence of a node, presented in the handle's local forward
    /// orientation.
    virtual string get_sequence(const handle_t& handle) const;
        
    /// Loop over all the handles to next/previous (right/left) nodes. Passes
    /// them to a callback which returns false to stop iterating and true to
    /// continue. Returns true if we finished and false if we stopped early.
    virtual bool follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
        
    // Copy over the template for nice calls
    using HandleGraph::follow_edges;
        
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual void for_each_handle(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
        
    // Copy over the template for nice calls
    using HandleGraph::for_each_handle;
        
    /// Return the number of nodes in the graph
    virtual size_t node_size() const;
        
    // We also have some extra functions
        
    /// Get the inward-facing start handle for this net graph. Useful when
    /// working with traversals.
    const handle_t& get_start() const;
        
    /// Get the outward-facing end handle for this net graph. Useful when
    /// working with traversals.
    const handle_t& get_end() const;
        
    /// Returns true if the given handle represents a meta-node for a child
    /// chain or unary snarl, and false if it is a normal node actually in
    /// the net graph snarl's contents.
    bool is_child(const handle_t& handle) const;
        
    /// Get the handle in the backing graph reading into the child chain or
    /// unary snarl in the orientation represented by this handle to a node
    /// representing a child chain or unary snarl.
    handle_t get_inward_backing_handle(const handle_t& child_handle) const;
        
protected:
    
    /// Make a NetGraph without filling in any of the child indexes.
    NetGraph(const Visit& start, const Visit& end, const HandleGraph* graph, bool use_internal_connectivity = false);
    
    /// Add a unary child snarl to the indexes.
    void add_unary_child(const Snarl* unary);
        
    /// Add a chain of one or more non-unary snarls to the index.
    void add_chain_child(const Chain& chain);
    
    // Save the backing graph
    const HandleGraph* graph;
        
    // And the start and end handles that bound the snarl we are working on.
    handle_t start;
    handle_t end;
        
    // Should we use the internal connectivity of chain nodes and unary
    // snarl nodes?
    bool use_internal_connectivity;
    
    // We keep the unary snarl boundaries, reading in with the unary snarl
    // contents to the right.
    unordered_set<handle_t> unary_boundaries;
        
    // We keep a map from handles that enter the ends of chains to the
    // reverse handles to their fronts. Whenever the backing graph tells us
    // to emit the one, we emit the other instead. This makes them look like
    // one big node.
    unordered_map<handle_t, handle_t> chain_end_rewrites;
        
    // We keep basically the reverse map, from chain start in chain forward
    // orientation to chain end in chain forward orientation. This lets us
    // find the edges off the far end of a chian.
    unordered_map<handle_t, handle_t> chain_ends_by_start;
        
    // Stores whether a chain or unary snarl, identified by the ID of its
    // start handle, is left-left, right-right, or left-right connected.
    unordered_map<id_t, tuple<bool, bool, bool>> connectivity;
        
};
    
/**
 * A structure to keep track of the tree relationships between Snarls and perform utility algorithms
 * on them
 */
class SnarlManager {
public:
        
    /// Construct a SnarlManager for the snarls returned by an iterator
    /// Also covers iterators of chains of snarls.
    template <typename SnarlIterator>
    SnarlManager(SnarlIterator begin, SnarlIterator end);
        
    /// Construct a SnarlManager for the snarls contained in an input stream
    SnarlManager(istream& in);
        
    /// Default constructor
    SnarlManager() = default;
        
    /// Destructor
    ~SnarlManager() = default;
        
    /// Cannot be copied because of all the internal pointer indexes
    SnarlManager(const SnarlManager& other) = delete;
    SnarlManager& operator=(const SnarlManager& other) = delete;
        
    /// Can be moved
    SnarlManager(SnarlManager&& other) = default;
    SnarlManager& operator=(SnarlManager&& other) = default;
        
    /// Returns a vector of pointers to the children of a Snarl.
    /// If given null, returns the top-level root snarls.
    const vector<const Snarl*>& children_of(const Snarl* snarl) const;
        
    /// Returns a pointer to the parent of a Snarl or nullptr if there is none
    const Snarl* parent_of(const Snarl* snarl) const;
        
    /// Returns the Snarl that a traversal points into at either the start
    /// or end, or nullptr if the traversal does not point into any Snarl.
    /// Note that Snarls store the end Visit pointing out of rather than
    /// into the Snarl, so they must be reversed to query it.
    const Snarl* into_which_snarl(int64_t id, bool reverse) const;
        
    /// Returns the Snarl that a Visit points into. If the Visit contains a
    /// Snarl rather than a node ID, returns a pointer the managed version
    /// of that snarl.
    const Snarl* into_which_snarl(const Visit& visit) const;
        
    /// Get the Chain that the given snarl participates in. Instead of
    /// asking this class to walk the chain for you, use ChainIterators on
    /// this chain.
    const Chain* chain_of(const Snarl* snarl) const;
        
    /// Return true if a Snarl is part of a nontrivial chain of more than
    /// one snarl.
    bool in_nontrivial_chain(const Snarl* here) const;
        
    /// Get all the snarls in all the chains under the given parent snarl.
    /// If the parent snarl is null, gives the top-level chains that connect and contain the top-level root snarls.
    /// Unary snarls and snarls in trivial chains will be presented as their own chains.
    /// Snarls are not necessarily oriented appropriately given their ordering in the chain.
    /// Useful for making a net graph.
    const deque<Chain>& chains_of(const Snarl* snarl) const;
        
    /// Get the net graph of the given Snarl's contents, using the given
    /// backing HandleGraph. If use_internal_connectivity is false, each
    /// chain and unary child snarl is treated as an ordinary node which is
    /// assumed to be only traversable from one side to the other.
    /// Otherwise, traversing the graph works like it would if you actually
    /// went through the internal graphs fo child snarls.
    NetGraph net_graph_of(const Snarl* snarl, const HandleGraph* graph, bool use_internal_connectivity = true) const;
        
    /// Returns true if snarl has no children and false otherwise
    bool is_leaf(const Snarl* snarl) const;
        
    /// Returns true if snarl has no parent and false otherwise
    bool is_root(const Snarl* snarl) const;
        
    /// Returns a reference to a vector with the roots of the Snarl trees
    const vector<const Snarl*>& top_level_snarls() const;
        
    /// Reverses the orientation of a snarl
    void flip(const Snarl* snarl);
        
    /// Add the given snarl to the SnarlManager as neither a root nor a
    /// child of any other snarl. The snarl must eventually either be added
    /// to a parent snarl or as a root snarl through add_chain() or it will
    /// not be visible.
    const Snarl* add_snarl(const Snarl& new_snarl);
        
    /// Add the given chain of snarls that have already been added with
    /// add_snarl(). Parents the chain to the given parent snarl, also added
    /// with add_snarl(). If the parent is null, makes the chain a root
    /// chain and all of its snarls root snarls. Note that the chains are
    /// allowed to be reallocated until the last child chain of a snarl is
    /// added.
    void add_chain(const Chain& new_chain, const Snarl* chain_parent);
        
    /// Returns the Nodes and Edges contained in this Snarl but not in any child Snarls (always includes the
    /// Nodes that form the boundaries of child Snarls, optionally includes this Snarl's own boundary Nodes)
    pair<unordered_set<Node*>, unordered_set<Edge*> > shallow_contents(const Snarl* snarl, VG& graph,
                                                                       bool include_boundary_nodes) const;
        
    /// Returns the Nodes and Edges contained in this Snarl, including those in child Snarls (optionally
    /// includes Snarl's own boundary Nodes)
    pair<unordered_set<Node*>, unordered_set<Edge*> > deep_contents(const Snarl* snarl, VG& graph,
                                                                    bool include_boundary_nodes) const;
        
    /// Look left from the given visit in the given graph and gets all the
    /// attached Visits to nodes or snarls.
    vector<Visit> visits_left(const Visit& visit, VG& graph, const Snarl* in_snarl) const;
        
    /// Look left from the given visit in the given graph and gets all the
    /// attached Visits to nodes or snarls.
    vector<Visit> visits_right(const Visit& visit, VG& graph, const Snarl* in_snarl) const;
        
    /// Returns a map from all Snarl boundaries to the Snarl they point into. Note that this means that
    /// end boundaries will be reversed.
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_boundary_index() const;
        
    /// Returns a map from all Snarl start boundaries to the Snarl they point into.
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_start_index() const;
        
    /// Returns a map from all Snarl end boundaries to the Snarl they point into. Note that this means that
    /// end boundaries will be reversed.
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_end_index() const;
        
    /// Execute a function on all top level sites
    void for_each_top_level_snarl(const function<void(const Snarl*)>& lambda) const;
        
    /// Execute a function on all sites in a preorder traversal
    void for_each_snarl_preorder(const function<void(const Snarl*)>& lambda) const;

    /// Execute a function on all sites (under a given root if specified) in a postorder traversal
    void for_each_snarl_postorder(const function<void(const Snarl*)>& lambda,
                                  const Snarl* root = NULL) const;
        
    /// Execute a function on all top level sites in parallel
    void for_each_top_level_snarl_parallel(const function<void(const Snarl*)>& lambda) const;
        
    /// Execute a function on all sites in parallel
    void for_each_snarl_parallel(const function<void(const Snarl*)>& lambda) const;
        
    /// Given a Snarl that we don't own (like from a Visit), find the
    /// pointer to the managed copy of that Snarl.
    const Snarl* manage(const Snarl& not_owned) const;
        
private:
    
    /// Define the key type
    using key_t = pair<pair<int64_t, bool>, pair<int64_t, bool>>;
        
    /// Master list of the snarls in the graph.
    /// Use a deque so pointers never get invalidated but we still have some locality.
    deque<Snarl> snarls;
        
    /// Roots of snarl trees
    vector<const Snarl*> roots;
    /// Chains of root-level snarls. Uses a deque so Chain* pointers don't get invalidated.
    deque<Chain> root_chains;
        
    /// Map of snarls to the child snarls they contain
    unordered_map<key_t, vector<const Snarl*>> children;
    /// Map of snarls to the child chains they contain
    /// Uses a deque so Chain* pointers don't get invalidated.
    unordered_map<key_t, deque<Chain>> child_chains;
    /// Map of snarls to their parent snarls
    unordered_map<key_t, const Snarl*> parent;
    /// Map of snarls to the chain each appears in
    unordered_map<key_t, const Chain*> parent_chain;
        
    /// Map of snarl keys to the pointer to the managed copy in the snarls vector.
    /// Is non-const so we can do flip nicely.
    unordered_map<key_t, Snarl*> self;
        
    /// Map of node traversals to the snarls they point into
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_into;
        
    /// Converts Snarl to the form used as keys in internal data structures
    inline key_t key_form(const Snarl* snarl) const;
        
    /// Builds tree indexes after Snarls have been added to the snarls vector
    void build_indexes();
        
    /// Actually compute chains for a set of already indexed snarls, which
    /// is important when chains were not provided. Returns the chains.
    deque<Chain> compute_chains(const vector<const Snarl*>& input_snarls);
        
    // Chain computation uses these pseudo-chain-traversal functions, which
    // walk around based on the snarl boundary index. This basically gets
    // you chains, except for structures that look like circular chains,
    // which are actually presented as circular chains and not linear ones.
    // They also let you walk into unary snarls.
        
    /// Get a Visit to the snarl coming after the given Visit to a snarl, or
    /// a Visit with no Snarl no next snarl exists. Accounts for snarls'
    /// orientations.
    Visit next_snarl(const Visit& here) const;
        
    /// Get a Visit to the snarl coming before the given Visit to a snarl,
    /// or a Visit with no Snarl no previous snarl exists. Accounts for
    /// snarls' orientations.
    Visit prev_snarl(const Visit& here) const;
        
    /// Get the Snarl, if any, that shares this Snarl's start node as either
    /// its start or its end. Does not count this snarl, even if this snarl
    /// is unary. Basic operation used to traverse a chain. Caller must
    /// account for snarls' orientations within a chain.
    const Snarl* snarl_sharing_start(const Snarl* here) const;
        
    /// Get the Snarl, if any, that shares this Snarl's end node as either
    /// its start or its end. Does not count this snarl, even if this snarl
    /// is unary. Basic operation used to traverse a chain. Caller must
    /// account for snarls' orientations within a chain.
    const Snarl* snarl_sharing_end(const Snarl* here) const;
};
    
/// Converts a Visit to a NodeTraversal. Throws an exception if the Visit is of a Snarl instead
/// of a Node
inline NodeTraversal to_node_traversal(const Visit& visit, const VG& graph);
    
/// Converts a Visit to a NodeTraversal in the opposite orientation. Throws an exception if the
/// Visit is of a Snarl instead of a Node
inline NodeTraversal to_rev_node_traversal(const Visit& visit, const VG& graph);
    
/// Converts a Visit to a node or snarl into a NodeSide for its left side.
inline NodeSide to_left_side(const Visit& visit);
    
/// Converts a Visit to a node or snarl into a NodeSide for its right side.
inline NodeSide to_right_side(const Visit& visit);
    
/// Converts a NodeTraversal to a Visit.
inline Visit to_visit(const NodeTraversal& node_traversal);
    
/// Converts a Mapping to a Visit. The mapping must represent a full node
/// match. If make_full_node_match is true, the mapping will automatically
/// be made a full node match during the conversion process.
inline Visit to_visit(const Mapping& mapping, bool make_full_node_match = false);
    
/// Make a Visit from a node ID and an orientation
inline Visit to_visit(id_t node_id, bool is_reverse);
    
/// Make a Visit from a snarl to traverse
inline Visit to_visit(const Snarl& snarl);
    
/// Get the reversed version of a visit
inline Visit reverse(const Visit& visit);
    
/// Converts a NodeTraversal to a Visit in the opposite orientation.
inline Visit to_rev_visit(const NodeTraversal& node_traversal);
    
/// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
/// of a Node. Uses a function to get node length.
inline Mapping to_mapping(const Visit& visit, std::function<size_t(id_t)> node_length);
    
/// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
/// of a Node. Uses a graph to get node length.
inline Mapping to_mapping(const Visit& visit, VG& vg);
    
/// Copies the boundary Visits from one Snarl into another
inline void transfer_boundary_info(const Snarl& from, Snarl& to);
    
// We need some Visit operators
    
/**
 * Two Visits are equal if they represent the same traversal of the same
 * Node or Snarl.
 */
bool operator==(const Visit& a, const Visit& b);
/**
 * Two Visits are unequal if they are not equal.
 */
bool operator!=(const Visit& a, const Visit& b);
/**
 * A Visit is less than another Visit if it represents a traversal of a
 * smaller node, or it represents a traversal of a smaller snarl, or it
 * represents a traversal of the same node or snarl forward instead of
 * backward.
 */
bool operator<(const Visit& a, const Visit& b);
    
/**
 * A Visit can be printed.
 */
ostream& operator<<(ostream& out, const Visit& visit);
    
// And some operators for SnarlTraversals
    
/**
 * Two SnarlTraversals are equal if their snarls are equal and they have the
 * same number of visits and all their visits are equal.
 */
bool operator==(const SnarlTraversal& a, const SnarlTraversal& b);
/**
 * Two SnarlTraversals are unequal if they are not equal.
 */
bool operator!=(const SnarlTraversal& a, const SnarlTraversal& b);
/**
 * A SnalTraversal is less than another if it is a traversal of a smaller
 * Snarl, or if its list of Visits has a smaller Visit first, or if its list
 * of Visits is shorter.
 */
bool operator<(const SnarlTraversal& a, const SnarlTraversal& b);
    
// And some operators for Snarls
    
/**
 * Two Snarls are equal if their types are equal and their bounding Visits
 * are equal and their parents are equal.
 */
bool operator==(const Snarl& a, const Snarl& b);
/**
 * Two Snarls are unequal if they are not equal.
 */
bool operator!=(const Snarl& a, const Snarl& b);
/**
 * A Snarl is less than another Snarl if its type is smaller, or its start
 * Visit is smaller, or its end Visit is smaller, or its parent is smaller.
 */
bool operator<(const Snarl& a, const Snarl& b);
    
/**
 * A Snarl can be printed.
 */
ostream& operator<<(ostream& out, const Snarl& snarl);
    
/****
 * Template and Inlines:
 ****/
    
template <typename SnarlIterator>
SnarlManager::SnarlManager(SnarlIterator begin, SnarlIterator end) {
    // add snarls to master list
    for (auto iter = begin; iter != end; iter++) {
        snarls.push_back(*iter);
    }
    // record the tree structure and build the other indexes
    build_indexes();
}
    
inline NodeTraversal to_node_traversal(const Visit& visit, VG& graph) {
    assert(visit.node_id());
    return NodeTraversal(graph.get_node(visit.node_id()), visit.backward());
}
    
inline NodeTraversal to_rev_node_traversal(const Visit& visit, VG& graph) {
    assert(visit.node_id());
    return NodeTraversal(graph.get_node(visit.node_id()), !visit.backward());
}
    
inline NodeSide to_left_side(const Visit& visit) {
    assert(visit.node_id() || (visit.snarl().start().node_id() && visit.snarl().end().node_id()));
    if (visit.node_id()) {
        // Just report the left side of this node
        return NodeSide(visit.node_id(), visit.backward());
    } else if (visit.backward()) {
        // This is a reverse visit to a snarl, so its left side is the right
        // side of the end visit of the snarl.
        assert(visit.snarl().end().node_id());
        return to_right_side(visit.snarl().end());
    } else {
        // This is a forward visit to a snarl, so its left side is the left
        // side of the start visit of the snarl.
        assert(visit.snarl().start().node_id());
        return to_left_side(visit.snarl().start());
    }
}
    
inline NodeSide to_right_side(const Visit& visit) {
    assert(visit.node_id() || (visit.snarl().start().node_id() && visit.snarl().end().node_id()));
    if (visit.node_id()) {
        // Just report the right side of this node
        return NodeSide(visit.node_id(), !visit.backward());
    } else if (visit.backward()) {
        // This is a reverse visit to a snarl, so its right side is the
        // left side of the start visit of the snarl.
        assert(visit.snarl().start().node_id());
        return to_left_side(visit.snarl().start());
    } else {
        // This is a forward visit to a snarl, so its right side is the
        // right side of the end visit of the snarl.
        assert(visit.snarl().end().node_id());
        return to_right_side(visit.snarl().end());
    }
}
    
inline Visit to_visit(const NodeTraversal& node_traversal) {
    Visit to_return;
    to_return.set_node_id(node_traversal.node->id());
    to_return.set_backward(node_traversal.backward);
    return to_return;
}
    
inline Visit to_visit(const Mapping& mapping, bool make_full_node_match) {
    if (!make_full_node_match) {
        // If we're not explicitly coercing the mapping to a full node match, make sure it already is one.
        assert(mapping_is_match(mapping));
        assert(mapping.position().offset() == 0);
    }
    Visit to_return;
    to_return.set_node_id(mapping.position().node_id());
    to_return.set_backward(mapping.position().is_reverse());
    return to_return;
}
    
inline Visit to_visit(id_t node_id, bool is_reverse) {
    Visit to_return;
    to_return.set_node_id(node_id);
    to_return.set_backward(is_reverse);
    return to_return;
}
    
inline Visit to_visit(const Snarl& snarl) {
    Visit to_return;
    // Only copy necessary fields
    *to_return.mutable_snarl()->mutable_start() = snarl.start();
    *to_return.mutable_snarl()->mutable_end() = snarl.end();
    return to_return;
}
    
inline Visit reverse(const Visit& visit) {
    // Copy the visit
    Visit to_return = visit;
    // And flip its orientation bit
    to_return.set_backward(!visit.backward());
    return to_return;
}
    
inline Visit to_rev_visit(const NodeTraversal& node_traversal) {
    Visit to_return;
    to_return.set_node_id(node_traversal.node->id());
    to_return.set_backward(!node_traversal.backward);
    return to_return;
}
    
inline Mapping to_mapping(const Visit& visit, std::function<size_t(id_t)> node_length) {
    // Can't have a Mapping to a snarl
    assert(visit.node_id());
    
    // Make a mapping to the right place
    Mapping mapping;
    mapping.mutable_position()->set_node_id(visit.node_id());
    mapping.mutable_position()->set_is_reverse(visit.backward());
        
    // Get the length of the node visited
    size_t length = node_length(visit.node_id());
        
    // Fill the Mapping in as a perfect match of that lenght
    Edit* match = mapping.add_edit();
    match->set_from_length(length);
    match->set_to_length(length);
        
    return mapping;
}
    
inline Mapping to_mapping(const Visit& visit, VG& graph) {
    return to_mapping(visit, [&](id_t id) {
            return graph.get_node(id)->sequence().size();
        });
}
    
inline void transfer_boundary_info(const Snarl& from, Snarl& to) {
    *to.mutable_start() = from.start();
    *to.mutable_end() = from.end();
}
    
}

// note: this hash funtion is not used internally because we want the internal indices to ignore any
// additional information so that the Snarls stored as references map to the same place
// as the original objects
namespace std {
/// hash function for Snarls
template<>
struct hash<const vg::Snarl> {
    size_t operator()(const vg::Snarl& snarl) const {
        auto hsh = hash<pair<pair<int64_t, bool>, pair<int64_t, bool> > >();
        return hsh(make_pair(make_pair(snarl.start().node_id(), snarl.start().backward()),
                             make_pair(snarl.end().node_id(), snarl.end().backward())));
    }
};
}

#endif /* snarls_hpp */
