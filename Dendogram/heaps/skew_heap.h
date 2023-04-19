#pragma once

#include "gbbs/gbbs.h"

namespace gbbs{
namespace skew_heap{
template <typename W>
struct skew_heap_node{
    W key;
    size_t val;
    skew_heap_node* left;
    skew_heap_node* right;
    
    skew_heap_node(W _key, size_t _val, skew_heap_node* _left = nullptr, skew_heap_node* _right = nullptr) 
    : key(_key), val(_val), left(_left), right(_right) {}
};

template <typename W>
skew_heap_node<W>* meld(skew_heap_node<W>* L, skew_heap_node<W>* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->key < R->key || (L->key == R->key && L->val <= R->val)){
            if (L->left == nullptr){
                L->left = R;
            } else{
                skew_heap_node<W>* temp = L->left;
                L->left = meld(L->right, R);
                L->right = temp;
            }
            return L;
        } else {
            return meld(R, L);
        }
    }
}

template <typename W>
struct skew_heap{
    skew_heap_node<W>* root;
    size_t size;

    skew_heap(skew_heap_node<W>* _root = nullptr)
      : root(_root){
        if (root){
            size = 1;
        } else{
            size = 0;
        }
      }
    
    std::pair<W, size_t> find_min(){
        return {root->key, root->val};
    }

    void delete_min(){
        skew_heap_node<W>* temp = root;
        root = meld(root->left, root->right);
        delete temp;
        size--;
    }

    bool is_empty(){
        return (size==0);
    }

    void insert(size_t val, W key){
        auto new_node = new skew_heap_node<W>(key, val);
        root = meld(root, new_node);
        size++;
    }

    void merge(skew_heap<W>* heap){
        root = meld(root, heap->root);
        size += heap->size;
        delete heap;
    }
};

}
}