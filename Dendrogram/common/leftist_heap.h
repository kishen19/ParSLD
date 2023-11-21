#pragma once
#include "gbbs/gbbs.h"

namespace gbbs{

namespace leftist_heap{

template <typename key_type>
struct node{
    key_type key;
    node* left;
    node* right;
    size_t rank;

    node(key_type _key) : key(_key), left(nullptr), right(nullptr) {
        rank = 1;
    }

    void init(key_type _key) {
        key = _key; 
        left = nullptr;
        right = nullptr;
        rank = 1;
    }

    void reassign_rank(){
        if (right){
            rank = right->rank + 1;
        } else{
            rank = 1;
        }
    }
    
    void swap_child(){
        if (left->rank > right->rank){
            node* temp = left;
            left = right;
            right = temp;
        }
    }
};

template <typename key_type>
node<key_type>* meld(node<key_type>* L, node<key_type>* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->key <= R->key){
            if (L->left == nullptr){
                L->left = R;
            } else{
                L->right = meld(L->right, R);
                L->swap_child();
                L->reassign_rank();
            }
            return L;
        } else {
            if (R->left == nullptr){
                R->left = L;
            } else{
                R->right = meld(R->right, L);
                R->swap_child();
                R->reassign_rank();
            }
            return R;
        }
    }
}

/////////// Leftist Heap

template <typename key_type>
struct heap{
    using node_allocator = parlay::type_allocator<node<key_type>>;
    node<key_type>* root;
    size_t size;

    heap() : root(nullptr), size(0){}
    
    inline key_type get_min(){
        return root->key;
    }

    inline void delete_min(){
        if (root){
            // node<key_type>* temp = root;
            root = meld(root->left, root->right);
            // node_allocator::destroy(temp);
            size--;
        } else{
            std::cout << "Empty Heap!" << std::endl;
            std::exit(1);
        }
    }

    inline bool is_empty(){
        return (size == 0);
    }

    inline void insert(key_type key){
        auto new_node = node_allocator::create(key);
        root = meld(root, new_node);
        size++;
    }

    inline void merge(heap<key_type>* H){
        root = meld(root, H->root);
        size += H->size;
        H->root = nullptr;
        H->size = 0;
    }

    template <class Seq>
    inline void init(const Seq& A){
        root = heapify_dc(A);
        size = A.size();
    }

    inline void init(node<key_type>* A, size_t n){
        root = heapify_dc(A, n);
        size = n;
    }

    template <class Seq>
    node<key_type>* heapify_dc(const Seq& A){
        size_t n = A.size();
        if (n == 0) {return nullptr;}
        else if (n <= 2048){
            // std::cout << "here, size= "<<n << std::endl;
            auto temp = node_allocator::create(A[0]);
            for (size_t i=1; i<n; i++){
                auto new_node = node_allocator::create(A[i]);
                temp = meld(temp, new_node);
            }
            return temp;
        } else {
            node<key_type> *heap1, *heap2;
            parlay::par_do(
                [&](){heap1 = heapify_dc(A.cut(0, n/2)); },
                [&](){heap2 = heapify_dc(A.cut(n/2, n)); }
            );
            return meld(heap1, heap2);
        }
    }

    node<key_type>* heapify_dc(node<key_type>* A, size_t n){
        if (n == 0) {return nullptr;}
        else if (n <= 2048){
            // std::cout << "here, size= "<<n << std::endl;
            auto temp = &A[0];
            for (size_t i=1; i<n; i++){
                temp = meld(temp, &A[i]);
            }
            return temp;
        } else {
            node<key_type> *heap1, *heap2;
            parlay::par_do(
                [&](){heap1 = heapify_dc(A, n/2); },
                [&](){heap2 = heapify_dc(A+n/2, n-n/2); }
            );
            return meld(heap1, heap2);
        }
    }
};

} // namespace leftist_heap
} // namespace gbbs