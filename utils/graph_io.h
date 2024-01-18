#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/graph_io.h"

namespace gbbs {

template <class Graph>
void write_sym_graph_binary(Graph& GA, char* outf, size_t n_batches = 8) {
  std::ofstream out(outf, std::ofstream::out | std::ios::binary);
  size_t n = GA.n;
  size_t m = GA.m;
  using W = typename Graph::weight_type;
  using edge_type = std::tuple<uintE, W>;

  // 1. Calculate total size
  auto offsets = parlay::sequence<uintT>(n + 1);
  std::cout << "# calculating size" << std::endl;
  parallel_for(0, n,
               [&](size_t i) { offsets[i] = GA.get_vertex(i).out_degree(); });
  offsets[n] = 0;
  size_t offset_scan = parlay::scan_inplace(make_slice(offsets));
  std::cout << "# offset_scan = " << offset_scan << " m = " << m << std::endl;
  assert(offset_scan == m);
  offsets[n] = m;

  long* sizes = gbbs::new_array_no_init<long>(3);
  std::cout << "Wrote: GA.n = " << GA.n << " sizeof edge = " << sizeof(edge_type) << std::endl;
  sizes[0] = GA.n;
  sizes[1] = GA.m;
  sizes[2] = sizeof(long) * 3 + sizeof(uintT) * (n + 1) + sizeof(edge_type) * m;
  out.write((char*)sizes, sizeof(long) * 3);   // write n, m and space used
  out.write((char*)offsets.begin(), sizeof(uintT) * (n + 1));   // write offsets

  // 2. Create binary format in-memory (batched)
  size_t block_size = parlay::num_blocks(n, n_batches);
  size_t edges_written = 0;
  for (size_t b = 0; b < n_batches; b++) {
    size_t start = b * block_size;
    size_t end = std::min(start + block_size, n);
    if (start >= end)
      break;
    std::cout << "# writing vertices " << start << " to " << end << std::endl;

    // create slab of graph and write out
    size_t start_offset = offsets[start];
    size_t end_offset = offsets[end];
    size_t n_edges = end_offset - start_offset;
    edge_type* edges = gbbs::new_array_no_init<edge_type>(n_edges);

    parallel_for(start, end, [&](size_t i) {
      size_t our_offset = offsets[i] - start_offset;
      auto map_f = [&](const uintE& u, const uintE& v, const W& wgh) {
        return wgh;
      };
      auto write_f = [&](const uintE& ngh, const uintT& offset, const W& val) {
        edges[offset] = std::make_tuple(ngh, val);
      };
      GA.get_vertex(i).out_neighbors().copy(our_offset, map_f, write_f);
    });

    size_t edge_space = sizeof(edge_type) * n_edges;
    out.write((char*)edges, edge_space);   // write edges

    std::cout << "# finished writing vertices " << start << " to " << end
              << std::endl;
    edges_written += n_edges;

    gbbs::free_array(edges, n_edges);
  }
  std::cout << "# Wrote " << edges_written << " edges in total, m = " << m
            << std::endl;
  assert(edges_written == m);
  out.close();

  // print_graph(GA);
}

}  // namespace gbbs

