#pragma once

namespace gbbs {

// 16 bytes / struct
template <typename weight_type>
struct RCtree_node {
  uintE parent;       // the parent of this node
  uintE edge_index;   // the edge stored at this node
  weight_type wgh;    // the weight of this edge.
  uintE round;   // the round this node is created
};

}   // namespace gbbs
