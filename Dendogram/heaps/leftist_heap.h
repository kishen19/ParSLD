#pragma once

#include "gbbs/gbbs.h"

namespace gbbs{

namespace leftist_heap{
template <typename W>
struct leftist_heap_node{
    W key;
    size_t val;
    leftist_heap_node* left;
    leftist_heap_node* right;
    size_t rank;

    void reassign_rank(){
        if (right){
            rank = right->rank + 1;
        } else{
            rank = 1;
        }
    }
    
    leftist_heap_node(W _key, size_t _val, leftist_heap_node* _left = nullptr, leftist_heap_node* _right = nullptr) 
    : key(_key), val(_val), left(_left), right(_right) {
        reassign_rank();
    }
    
    void swap_child(){
        if (left->rank > right->rank){
            leftist_heap_node* temp = left;
            left = right;
            right = temp;
        }
    }
};

template <typename W>
leftist_heap_node<W>* meld(leftist_heap_node<W>* L, leftist_heap_node<W>* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->key < R->key || (L->key == R->key && L->val <= R->val)){
            if (L->left == nullptr){
                L->left = R;
            } else{
                L->right = meld(L->right, R);
                L->swap_child();
                L->reassign_rank();
            }
            return L;
        } else {
            return meld(R, L);
        }
    }
}

template <typename W>
struct leftist_heap{
    leftist_heap_node<W>* root;
    size_t size;

    leftist_heap(leftist_heap_node<W>* _root = nullptr)
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
        if (root){
            leftist_heap_node<W>* temp = root;
            root = meld(root->left, root->right);
            delete temp;
            size--;
        }
    }

    bool is_empty(){
        return (size==0);
    }

    void insert(size_t val, W key){
        auto new_node = new leftist_heap_node<W>(key, val);
        root = meld(root, new_node);
        size++;
    }

    void merge(leftist_heap<W>* heap){
        root = meld(root, heap->root);
        size += heap->size;
        heap->root = nullptr;
        heap->size = 0;
    }
};

}
}