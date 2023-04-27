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
};

}
}