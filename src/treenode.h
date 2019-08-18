#ifndef TREENODE_HH
#define TREENODE_HH
#define UIDLEN 8

#include <iostream>
#include <stdio.h>
#include <exception>
#include <vector>

class IndexOutOfBoundsException: public std::exception
{
    virtual const char* what() const throw()
    {
        return "Index out of bounds.";
    };
};

/**
 * @brief treenode class
 * 
 */
class treenode {

private:

    /**
     * @brief Unique node id
     * 
     */
    std::string uid;

    /**
     * @brief The node's children
     * 
     */
    std::vector<treenode *> children;

    /**
     * @brief The node's parent node
     * 
     */
    treenode *parent;

    /**
     * @brief Generates random uid
     * 
     * @return std::string& 
     */
    const static std::string& generate_random_uid();

public:
    typedef std::vector<treenode *> tree_t;
    typedef tree_t::const_iterator tree_const_it;

    treenode();

    /**
     * @brief Construct a new treenode object
     * 
     * @param &uid Unique id
     * @param *parent Parent node
     * @param children The node's children
     */
    treenode( const std::string& uid, treenode *parent, std::vector<treenode *>& children );

    ~treenode();

    /**
     * @brief Adds a node to the first position
     * 
     * @param *child The node to add
     */
    void insert_child( treenode *child );

    /**
     * @brief Adds a child at the given position
     * 
     * @param *child The node to add
     * @param &position The desired position
     */
    void add_child_at( treenode *child, const int& position );

    /**
     * @brief Appends child as the last child
     * 
     * @param *child The node to append
     * @return int& The child's position
     */
    const int& append_child( treenode *child );

    /**
     * @brief Removes a child node given its instance
     * 
     * @param *child Child to be removed
     * @return int& The position where the removed node was
     */
    const int& erase_child( treenode *child );

    /**
     * @brief Removes a node at the given position
     * 
     * @param &position Position where the node should be
     */
    void pop_child_at( const uint& position );

    /**
     * @brief Removes the last child
     * 
     * @return int& Position where the removed node was
     */
    const int& pop_child();

    /**
     * @brief Counts this node's children
     * 
     * @return const int& The children count
     */
    const int& children_count() const;

    const int& node_count() const;

    const bool& is_leaf() const;

    const bool& is_root() const;

    /**
     * Getters and setters
     */

    treenode* get_parent();

    void set_parent( treenode *parent );

    // TODO
    tree_const_it children_begin();

    // TODO
    tree_const_it children_at( const uint& position );

    // TODO
    tree_const_it children_end();

};

#endif