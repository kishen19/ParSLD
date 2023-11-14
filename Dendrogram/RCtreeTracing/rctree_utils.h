#pragma once

namespace gbbs {

size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

// ~32 bytes / struct
template <typename weight_type>
struct RCtree_node {
  size_t parent;       // the parent of this node
  size_t edge_index;   // the edge stored at this node
  size_t alt =
     SIZE_T_MAX;   // used in the compress case. (default = SIZE_T_MAX)
  uintE round = UINT_E_MAX;   // the round this node is created
  weight_type wgh;            // the weight of this edge.
};

}   // namespace gbbs
