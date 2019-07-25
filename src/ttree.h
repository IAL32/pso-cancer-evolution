#ifndef TTREE_H
#define TTREE_H

#include "tree.h"


/**
 * @brief Tumoral tree
 * 
 */
class ttree {
/**
 * @brief Derivating everything from the tree class
 * 
 */
friend class tree;

private:
    uint mutation_id;
    std::string name;
    bool loss;

public:

    // random tree operations

    static int add_random_back_mutation( tree& tree, uint k );
    
    static int delete_random_back_mutation( tree& tree );

    static int switch_random_nodes( tree& tree );

    static int 

};

#endif