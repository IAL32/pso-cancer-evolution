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

    /**
     * @brief Unique node id
     * 
     */
    uint uid;

    /**
     * @brief The node's children
     * 
     */
    std::vector<tree> children;

    /**
     * @brief The node's parent node
     * 
     */
    tree& parent;

    /**
     * @brief The node's first child
     * 
     */
    tree& first_child;

    /**
     * @brief The node's last child
     * 
     */
    tree& last_child;

public:

    tree();

    /**
     * @brief Construct a new tree object
     * 
     * @param uid Unique id
     * @param parent Parent node
     * @param children The node's children
     */
    tree( uint uid, tree& parent, std::vector<tree> children );

    ~tree();

    /**
     * @brief Adds a node to the first position
     * 
     * @param child The node to add
     */
    void add_child( tree& child );

    /**
     * @brief Adds a child at the given position
     * 
     * @param child The node to add
     * @param position The desired position
     */
    void add_child_at( tree& child, uint position );

    /**
     * @brief Appends child as the last child
     * 
     * @param child The node to append
     * @return int The child's position
     */
    int append_child( tree& child );

    /**
     * @brief Removes a child node given its instance
     * 
     * @param child Child to be removed
     * @return int The position where the removed node was
     */
    int remove_child( tree& child );

    /**
     * @brief Removes a node at the given position
     * 
     * @param position Position where the node should be
     */
    void remove_child_at( uint position );

    /**
     * @brief Removes the last child
     * 
     * @return int Position where the removed node was
     */
    int remove_last_child();

};

#endif