#include "position.hpp"

namespace vg {

pos_t make_pos_t(const Position& pos) {
    return make_tuple(pos.node_id(), pos.is_reverse(), pos.offset());
}

pos_t make_pos_t(id_t id, bool is_rev, off_t off) {
    return make_tuple(id, is_rev, off);
}

Position make_position(const pos_t& pos) {
    Position p;
    p.set_node_id(id(pos));
    p.set_is_reverse(is_rev(pos));
    p.set_offset(offset(pos));
    return p;
}

Position make_position(id_t id, bool is_rev, off_t off) {
    Position p;
    p.set_node_id(id);
    p.set_is_reverse(is_rev);
    p.set_offset(off);
    return p;
}

bool is_empty(const pos_t& pos) {
    return id(pos) == 0;
}

id_t id(const pos_t& pos) {
    return get<0>(pos);
}

bool is_rev(const pos_t& pos) {
    return get<1>(pos);
}

off_t offset(const pos_t& pos) {
    return get<2>(pos);
}

id_t& get_id(pos_t& pos) {
    return get<0>(pos);
}

bool& get_is_rev(pos_t& pos) {
    return get<1>(pos);
}

off_t& get_offset(pos_t& pos) {
    return get<2>(pos);
}

pos_t reverse(const pos_t& pos, size_t node_length) {
    pos_t rev = pos;
    // swap the offset onto the other strand
    get_offset(rev) = node_length - offset(rev);
    // invert the position
    get_is_rev(rev) = !is_rev(rev);
    return rev;
}

Position reverse(const Position& pos, size_t node_length) {
    auto p = pos;
    p.set_offset(node_length - pos.offset());
    p.set_is_reverse(!pos.is_reverse());
    return p;
}

ostream& operator<<(ostream& out, const pos_t& pos) {
    return out << id(pos) << (is_rev(pos) ? "-" : "+") << offset(pos);
}

size_t xg_cached_node_length(id_t id, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache) {
    //cerr << "Looking for position " << pos << endl;
    pair<Node, bool> cached = node_cache.retrieve(id);
    if(!cached.second) {
        //cerr << "Not in the cache" << endl;
        // If it's not in the cache, put it in
        cached.first = xgidx->node(id);
        node_cache.put(id, cached.first);
    }
    Node& node = cached.first;
    return node.sequence().size();
}

char xg_cached_pos_char(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache) {
    //cerr << "Looking for position " << pos << endl;
    pair<Node, bool> cached = node_cache.retrieve(id(pos));
    if(!cached.second) {
        //cerr << "Not in the cache" << endl;
        // If it's not in the cache, put it in
        cached.first = xgidx->node(id(pos));
        node_cache.put(id(pos), cached.first);
    }
    Node& node = cached.first;
    if (is_rev(pos)) {
        /*
        cerr << "reversed... " << endl;
        cerr << "rev pos " << offset(reverse(pos, node.sequence().size())) << endl;
        cerr << "seq is " << node.sequence() << " and got " <<
            reverse_complement(node.sequence()[offset(reverse(pos, node.sequence().size()))-1]) << endl;
        */
        return reverse_complement(node.sequence()[offset(reverse(pos, node.sequence().size()))-1]);
    } else {
        /*
        cerr << "forward... " << endl;
        cerr << "seq is " << node.sequence() << " and got " << node.sequence().at(offset(pos)) << endl;
        */
        return node.sequence().at(offset(pos));
    }
}

map<pos_t, char> xg_cached_next_pos_chars(pos_t pos, xg::XG* xgidx, LRUCache<id_t, Node>& node_cache) {

    map<pos_t, char> nexts;
    // See if the node is cached (did we just visit it?)
    pair<Node, bool> cached = node_cache.retrieve(id(pos));
    if(!cached.second) {
        // If it's not in the cache, put it in
        cached.first = xgidx->node(id(pos));
        node_cache.put(id(pos), cached.first);
    }
    Node& node = cached.first;
    // if we are still in the node, return the next position and character
    if (offset(pos) < node.sequence().size()-1) {
        ++get_offset(pos);
        nexts[pos] = xg_cached_pos_char(pos, xgidx, node_cache);
    } else {

        auto is_inverting = [](const Edge& e) {
            return !(e.from_start() == e.to_end())
            && (e.from_start() || e.to_end());
        };

        // look at the next positions we could reach

        if (!is_rev(pos)) {
            // we are on the forward strand, the next things from this node come off the end
            for (auto& edge : xgidx->edges_on_end(id(pos))) {
                id_t nid = (edge.from() == id(pos) ?
                            edge.to()
                            : edge.from());
                pos_t p = make_pos_t(nid, is_inverting(edge), 0);
                nexts[p] = xg_cached_pos_char(p, xgidx, node_cache);
            }
        } else {
            // we are on the reverse strand, the next things from this node come off the start
            for (auto& edge : xgidx->edges_on_start(id(pos))) {
                id_t nid = (edge.to() == id(pos) ?
                            edge.from()
                            : edge.to());
                pos_t p = make_pos_t(nid, !is_inverting(edge), 0);
                nexts[p] = xg_cached_pos_char(p, xgidx, node_cache);
            }
        }
    }
    return nexts;
}

}
