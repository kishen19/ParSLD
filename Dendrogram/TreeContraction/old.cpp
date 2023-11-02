template <class Graph>
void DendrogramContract(Graph& GA){
	using W = typename Graph::weight_type;
	// using edge = std::tuple(uintE, uintE);
	// using edge_ind = std::tuple(uintE, uintE, size_t);
	// using heap = leftist_heap::heap<kv>;
	// using heap = skew_heap::heap<kv>;
	using heap = pairing_heap::heap<kv>;
	size_t n = GA.n;
	size_t m = GA.m;

    // auto heaps = gbbs::new_array<heap>(n);
	auto uf = simple_union_find::SimpleUnionAsyncStruct(n);
	auto deg = sequence<size_t>::from_function(n, [&](uintE i){
		return GA.get_vertex(i).out_degree();
	});
	auto deg_old = sequence<size_t>::uninitialized(n);
	for(uintE i=0; i<n; i++){
		std::cout << deg[i] << " ";
	}
	std::cout << std::endl;

	parlay::random_generator gen;
	std::uniform_int_distribution<> dis(0, 1);
	auto tosses = sequence<bool>(n,0);

	size_t rem = m;
	while (rem){
		parlay::copy(deg, deg_old);
		// Step 1: Rake
		parallel_for(0, n, [&](uintE i){
			if (deg_old[i]>1){
				auto filter_f = [&](const uintE& src, const uintE& nbhr, const W& wgh){
					return deg_old[nbhr] == 1;
				};
				auto out_f = [&](size_t j, const std::tuple<uintE, W>& nw){
					auto v = std::get<0>(nw);
					uf.unite(i,v);
					deg[v] -=1;
				};
				auto num = GA.get_vertex(i).out_neighbors().count(filter_f);
				auto tmp = sequence<std::tuple<uintE, W>>::uninitialized(num);
				GA.get_vertex(i).out_neighbors().filter(filter_f, out_f, tmp.begin());
				rem -= (2*num);
				deg[i] -= num;
			}
		});
		std::cout << rem << std::endl;
		parlay::copy(deg, deg_old);

		// Step 2: Compress
		// auto ind = parlay::iota(n);
		// auto deg2 = parlay::filter(ind, [&](size_t u){
		// 	return deg[u] == 2;
		// });
		// parallel_for(0, deg2.size(),[&](size_t i){
		// 	auto r = gen[i];
		// 	tosses[deg2[i]] =  dis(r);
		// });
		// gen += n;
		// parallel_for(0, deg2.size(), [&](size_t i){
		// 	auto u = deg2[i];
			
		// });


		// for(auto i=0; i < mis.size(); i++){
		// 	std::cout << mis[i] << " ";
		// }
		// std::cout << std::endl;
	}
}