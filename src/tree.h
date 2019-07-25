#ifndef TREE_HH
#define TREE_HH

#include <iostream>
#include <stdio.h>
#include <vector>

class tree;

/**
 * @brief Tree class
 * 
 */
class tree {

private:
    uint uid;
    std::vector<tree> children;
    tree& parent;
    tree& first_child;
    tree& last_child;

public:

    tree();
    tree( uint uid, tree& parent, tree& first_child, tree& last_child );

};

#endif