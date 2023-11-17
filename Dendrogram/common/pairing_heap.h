#pragma once

#include "gbbs/gbbs.h"
#include "binary_heap.h"

namespace gbbs{
namespace pairing_heap{

template <class Node>
Node* meld(Node* L, Node* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->get_key() <= R->get_key()){
            R->sibling = L->child;
            L->child = R;
            return L;
        } else {
            L->sibling = R->child;
            R->child = L;
            return R;
        }
    }
}

template <class Node>
Node* two_pass_merge(Node* H){
    if (H == nullptr || H->sibling == nullptr){
        return H;
    } else{
        Node* temp1 = H->sibling;
        Node* temp2 = H->sibling->sibling;
        H->sibling = nullptr;
        temp1->sibling = nullptr;
        return meld(meld(H, temp1), two_pass_merge(temp2));
    }
}

/////////// Pairing Heap

template <typename key_type>
struct node{
    key_type key;
    node* child;
    node* sibling;
    
    node(key_type _key) : key(_key), child(nullptr), sibling(nullptr) {}

    void init(key_type _key) {
        key = _key; 
        child = nullptr; 
        sibling = nullptr;
    }

    inline key_type get_key(){
        return key;
    }
};

template <typename key_type>
struct heap{
    using node_allocator = parlay::type_allocator<node<key_type>>;
    node<key_type>* root;
    size_t size;

    heap() : root(nullptr), size(0){}
    
    inline key_type get_min(){
        if (root){
            return root->key;
        } else{
            std::cout << "Empty Heap!" << std::endl;
            std::exit(1);
        }
    }

    inline void delete_min(){
        if (root){
            node<key_type>* temp = root;
            root = two_pass_merge(root->child);
            node_allocator::destroy(temp);
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
    size_t filter(key_type key, Seq& out){ //sequential
        size_t k=0;
        while (root && root->get_key() <= key){
            auto kv = get_min();
            delete_min();
            out[k] = kv;
            k++;
        }
        size -= k;
        return k;
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

/////////// Block Pairing Heap

template <typename key_type>
struct block_node{
    using node_allocator = parlay::type_allocator<binary_heap::heap<key_type>>;
    binary_heap::heap<key_type>* block;
    block_node* child;
    block_node* sibling;

    block_node(key_type* keys, size_t num) : child(nullptr), sibling(nullptr) {
        block = node_allocator::create(keys, num);
    }

    ~block_node(){
        node_allocator::destroy(block);
    }

    inline key_type get_key(){
        return block->get_min();
    }
};

template <typename key_type>
struct block_heap{
    using node_allocator = parlay::type_allocator<block_node<key_type>>;
    block_node<key_type>* root;

    inline void init(key_type* block, size_t num) {
        if (num > 0){
            root = node_allocator::create(block, num);
        } else{
            root = nullptr;
        }
    }
    
    inline key_type get_min(){
        if (root){
            return root->block->get_min();
        } else{
            std::cout << "Empty Heap" << std::endl;
            std::exit(1);
        }
    }

    void delete_min(){
        if (root){
            root->block->delete_min();

            if (root->block->size > 0){
                root->sibling = root->child;
                root->child = nullptr;
            } else{
                block_node<key_type>* temp = root;
                root = root->child;
                node_allocator::destroy(temp);
            }
            root = two_pass_merge(root);
        } else{
            std::cout << "Empty Heap" << std::endl;
            std::exit(1);
        }
    }

    inline bool is_empty(){
        return (root == nullptr);
    }

    inline void merge(block_heap<key_type>* H){
        root = meld(root, H->root);
        H->root = nullptr;
    }
};

} // namespace pairing_heap
} // namespace gbbs