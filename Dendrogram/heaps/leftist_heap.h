#pragma once

#include "gbbs/gbbs.h"

namespace gbbs{

namespace leftist_heap{

template <typename key_type>
struct leftist_heap_node{
    key_type key;
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
    
    leftist_heap_node(key_type _key, leftist_heap_node* _left = nullptr, leftist_heap_node* _right = nullptr) 
    : key(_key), left(_left), right(_right) {
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

template <typename key_type>
leftist_heap_node<key_type>* meld(leftist_heap_node<key_type>* L, leftist_heap_node<key_type>* R){
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

template <typename key_type>
struct leftist_heap{
    leftist_heap_node<key_type>* root;

    leftist_heap(leftist_heap_node<key_type>* _root = nullptr)
      : root(_root){}
    
    key_type find_min(){
        return root->key;
    }

    void delete_min(){
        if (root){
            leftist_heap_node<key_type>* temp = root;
            root = meld(root->left, root->right);
            delete temp;
        }
    }

    bool is_empty(){
        return (root == nullptr);
    }

    void insert(key_type key){
        leftist_heap_node<key_type>* new_node = new leftist_heap_node<key_type>(key);
        root = meld(root, new_node);
    }

    void merge(leftist_heap<key_type>* heap){
        root = meld(root, heap->root);
        heap->root = nullptr;
    }
};

}
}