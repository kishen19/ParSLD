#pragma once

#include "gbbs/gbbs.h"

namespace gbbs{

namespace pairing_heap{

template <typename key_type>
struct pairing_heap_node{
    key_type key;
    pairing_heap_node* child;
    pairing_heap_node* sibling;
    
    pairing_heap_node(key_type _key, pairing_heap_node* _child = nullptr, pairing_heap_node* _sibling = nullptr) 
    : key(_key), child(_child), sibling(_sibling) {}
};

template <typename key_type>
pairing_heap_node<key_type>* meld(pairing_heap_node<key_type>* L, pairing_heap_node<key_type>* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->key <= R->key){
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

template <typename key_type>
pairing_heap_node<key_type>* two_pass_merge(pairing_heap_node<key_type>* heap){
    if (heap == nullptr || heap->sibling == nullptr){
        return heap;
    } else{
        pairing_heap_node<key_type>* temp1 = heap->sibling;
        pairing_heap_node<key_type>* temp2 = heap->sibling->sibling;
        heap->sibling = nullptr;
        temp1->sibling = nullptr;
        return meld(meld(heap, temp1), two_pass_merge(temp2));
    }
}

template <typename key_type>
struct pairing_heap{
    pairing_heap_node<key_type>* root;

    pairing_heap(pairing_heap_node<key_type>* _root = nullptr)
      : root(_root){}
    
    key_type find_min(){
        return root->key;
    }

    void delete_min(){
        if (root){
            pairing_heap_node<key_type>* temp = root;
            root = two_pass_merge(root->child);
            delete temp;
        }
    }

    bool is_empty(){
        return (root == nullptr);
    }

    void insert(key_type key){
        pairing_heap_node<key_type>* new_node = new pairing_heap_node<key_type>(key);
        root = meld(root, new_node);
    }

    void merge(pairing_heap<key_type>* heap){
        root = meld(root, heap->root);
        heap->root = nullptr;
    }

    template <class Seq>
    pairing_heap_node<key_type>* heapify_dc(const Seq& A){
        auto n = A.size();
        if (n == 0) {return nullptr;}
        else if (n <= 256){
            auto temp = new pairing_heap_node<key_type>(A[0]);
            for (size_t i=1; i<n; i++){
                auto new_node = new pairing_heap_node<key_type>(A[i]);
                temp = meld(temp, new_node);
            }
            return temp;
        } else {
            pairing_heap_node<key_type> *heap1, *heap2;
            parlay::par_do(
                [&](){heap1 = heapify_dc(A.cut(0, n/2));},
                [&](){heap2 = heapify_dc(A.cut(n/2, n));}
            );
            return meld(heap1, heap2);
        }
    }

    template <class Seq>
    pairing_heap_node<key_type>* heapify_linear(const Seq& A){ // assuming sorted by (weights,ind)
        auto n = A.size();
        if (n == 0) {return nullptr;}
        else{
            timer t;
            t.start();
            auto nodes = gbbs::sequence<pairing_heap_node<key_type>*>::uninitialized(n);
            parallel_for(0, n, [&](size_t i){
                nodes[i] = new pairing_heap_node<key_type>(A[i]);
            });
            parallel_for(0, n, [&](size_t i){
                if (2*i+1 < n){
                    nodes[i]->child = nodes[2*i+1];
                }
                if (2*i+2 < n){
                    nodes[i]->child->sibling = nodes[2*i+2];
                }
            });
            return nodes[0];
        }
    }

    template <class Seq>
    void create_heap(const Seq& A){
        root = heapify_linear(A);
        // root = heapify_dc(A);
    }
};

}
}