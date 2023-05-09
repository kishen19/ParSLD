#pragma once

#include "gbbs/gbbs.h"

namespace gbbs{

namespace skew_heap{

template <typename key_type>
struct skew_heap_node{
    key_type key;
    skew_heap_node* left;
    skew_heap_node* right;
    
    skew_heap_node(key_type _key, skew_heap_node* _left = nullptr, skew_heap_node* _right = nullptr) 
    : key(_key), left(_left), right(_right) {}
};

template <typename key_type>
skew_heap_node<key_type>* meld(skew_heap_node<key_type>* L, skew_heap_node<key_type>* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->key <= R->key){
            if (L->left == nullptr){
                L->left = R;
            } else{
                skew_heap_node<key_type>* temp = L->left;
                L->left = meld(L->right, R);
                L->right = temp;
            }
            return L;
        } else {
            if (R->left == nullptr){
                R->left = L;
            } else{
                skew_heap_node<key_type>* temp = R->left;
                R->left = meld(R->right, L);
                R->right = temp;
            }
            return R;
        }
    }
}

template <typename key_type>
struct skew_heap{
    skew_heap_node<key_type>* root;

    skew_heap(skew_heap_node<key_type>* _root = nullptr)
      : root(_root){}
    
    key_type find_min(){
        return root->key;
    }

    void delete_min(){
        if (root){
            skew_heap_node<key_type>* temp = root;
            root = meld(root->left, root->right);
            delete temp;
        }
    }

    bool is_empty(){
        return (root == nullptr);
    }

    void insert(key_type key){
        skew_heap_node<key_type>* new_node = new skew_heap_node<key_type>(key);
        root = meld(root, new_node);
    }

    void merge(skew_heap<key_type>* heap){
        root = meld(root, heap->root);
        heap->root = nullptr;
    }

    template <class Seq>
    skew_heap_node<key_type>* heapify(const Seq& A){ // assuming sorted by (weights,ind)
        auto n = A.size();
        if (n == 0) {return nullptr;}
        else{
            auto nodes = gbbs::sequence<skew_heap_node<key_type>*>::uninitialized(n);
            parallel_for(0, n, 256, [&](size_t i){
                nodes[i] = new skew_heap_node<key_type>(A[i]);
            });
            parallel_for(0, n, 256, [&](size_t i){
                if (2*i+1 < n){
                    nodes[i]->left = nodes[2*i+1];
                }
                if (2*i+2 < n){
                    nodes[i]->right = nodes[2*i+2];
                }
            });
            return nodes[0];
        }
    }

    // template <class Seq>
    // skew_heap_node<key_type>* heapify(const Seq& A){
    //     auto n = A.size();
    //     if (n == 0) {return nullptr;}
    //     else if (n <= 256){
    //         auto temp = new skew_heap_node<key_type>(A[0]);
    //         for (size_t i=1; i<n; i++){
    //             auto new_node = new skew_heap_node<key_type>(A[i]);
    //             temp = meld(temp, new_node);
    //         }
    //         return temp;
    //     } else {
    //         skew_heap_node<key_type> *heap1, *heap2;
    //         parlay::par_do(
    //             [&](){heap1 = heapify(A.cut(0, n/2));},
    //             [&](){heap2 = heapify(A.cut(n/2, n));}
    //         );
    //         return meld(heap1, heap2);
    //     }
    // }

    template <class Seq>
    void create_heap(const Seq& A){
        root = heapify(A);
    }
};

}
}