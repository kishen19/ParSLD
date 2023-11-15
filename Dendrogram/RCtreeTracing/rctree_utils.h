#pragma once

namespace gbbs {

// 20 bytes / struct
template <typename weight_type>
struct RCtree_node {
  uintE parent;       // the parent of this node
  uintE edge_index;   // the edge stored at this node
  uintE alt;   // used in the compress case. (default = SIZE_T_MAX)
  uintE round;   // the round this node is created
  weight_type wgh;            // the weight of this edge.
};

}   // namespace gbbs
