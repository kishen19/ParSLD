#pragma once

namespace gbbs{
namespace binary_heap{

template <typename key_type>
struct heap{
    key_type* nodes;
    size_t size;

    heap(key_type* _nodes, size_t _size): nodes(_nodes), size(_size){}

    void heapify(size_t i){
        size_t l = 2*i+1;
        size_t r = 2*i+2;
        size_t ind = i;
        if (l < size && nodes[l] < nodes[i]){
            ind = l;
        }
        if (r < size && nodes[r] < nodes[ind]){
            ind = r;
        }
        if (ind != i){
            key_type temp = nodes[i];
            nodes[i] = nodes[ind];
            nodes[ind]= temp;
            heapify(ind);
        }
    }

    inline key_type get_min(){
        if (size > 0){
            return nodes[0];
        } else{
            std::cout << "Empty Heap!" << std::endl;
            std::exit(1);
        }
    }

    inline void delete_min(){
        if(size == 1){
            size--;
        } else if(size > 1){
            nodes[0] = nodes[size-1];
            size--;
            heapify(0);
        } else{
            std::cout << "Empty Heap!" << std::endl;
            std::exit(1);
        }
    }
};

} // namespace binary_heap
} // namespace gbbs