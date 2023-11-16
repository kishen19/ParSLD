#pragma once

namespace gbbs {

// 20 bytes / struct
template <typename weight_type>
struct RCtree_node {
  uintE parent;       // the parent of this node
  uintE edge_index;   // the edge stored at this node
  uintE alt;   // used in the compress case. (default = UINT_E_MAX)
  uintE round;   // the round this node is created
  weight_type wgh;            // the weight of this edge.
  int is_ready;  // takes the value {0, 1, 2}
  // - denotes how many neighbors have lower priority than us.
  // - if is_ready = 2, and degree = 2, then we can compress
};

}   // namespace gbbs
