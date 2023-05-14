#pragma once

namespace gbbs{

namespace skew_heap{

template <class Node>
Node* meld(Node* L, Node* R){
    if (L == nullptr){
        return R;
    } else if (R == nullptr){
        return L;
    } else{
        if (L->get_key() <= R->get_key()){
            if (L->left == nullptr){
                L->left = R;
            } else{
                Node* temp = L->left;
                L->left = meld(L->right, R);
                L->right = temp;
            }
            return L;
        } else {
            if (R->left == nullptr){
                R->left = L;
            } else{
                Node* temp = R->left;
                R->left = meld(R->right, L);
                R->right = temp;
            }
            return R;
        }
    }
}

/////////// Skew Heap

template <typename key_type>
struct node{
    key_type key;
    node* left;
    node* right;
    
    node(key_type _key, node* _left = nullptr, node* _right = nullptr) 
    : key(_key), left(_left), right(_right) {}

    inline key_type get_key(){
        return key;
    }
};

template <typename key_type>
struct heap{
    using node_allocator = parlay::type_allocator<node<key_type>>;
    node<key_type>* root;

    inline void init(key_type* A, size_t num){
        root = heapify_linear(A, num);
        // root = heapify_dc(A, num);
    }
    
    inline key_type get_min(){
        return root->key;
    }

    inline void delete_min(){
        if (root){
            node<key_type>* temp = root;
            root = meld(root->left, root->right);
            node_allocator::destroy(temp);
        } else{
            std::cout << "Empty Heap!" << std::endl;
            std::exit(1);
        }
    }

    bool is_empty(){
        return (root == nullptr);
    }

    inline void insert(key_type key){
        auto new_node = node_allocator::create(key);
        root = meld(root, new_node);
    }

    inline void merge(heap<key_type>* H){
        root = meld(root, H->root);
        H->root = nullptr;
    }

    node<key_type>* heapify_linear(key_type* A, size_t n){ // assuming sorted by (weights,ind)
        if (n == 0) {return nullptr;}
        else{
            auto nodes = gbbs::sequence<node<key_type>*>::uninitialized(n);
            parallel_for(0, n, 2048, [&](size_t i){
                nodes[i] = node_allocator::create(A[i]);
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

    node<key_type>* heapify_dc(key_type* A, size_t n){
        if (n == 0) {return nullptr;}
        else if (n <= 2048){
            auto temp = node_allocator::create(A[0]);
            for (size_t i=1; i<n; i++){
                auto new_node = node_allocator::create(A[i]);
                temp = meld(temp, new_node);
            }
            return temp;
        } else {
            node<key_type> *heap1, *heap2;
            parlay::par_do(
                [&](){heap1 = heapify_dc(A, n/2);},
                [&](){heap2 = heapify_dc(A + n/2, n - n/2);}
            );
            return meld(heap1, heap2);
        }
    }
};

} // namespace skew_heap
} // namespace gbbs