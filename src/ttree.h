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

    ttree( uint uid, tree& parent, tree& first_child, tree& last_child, std::string mutation_name, int mutation_id );

    // random tree operations

    int random_add_back_mutation( uint k );
    
    int random_delete_back_mutation();

    int random_delete_nodes();

    int random_prune_regraft();

};

#endif