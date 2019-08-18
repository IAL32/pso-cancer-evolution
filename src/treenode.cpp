#include "treenode.h"
#include <random>
#include <algorithm>

// https://stackoverflow.com/a/12468109
const std::string& treenode::generate_random_uid()
{
    auto randchar = []() -> char
    {
        const char charset[] = "0123456789" "ABCDEF";
        const uint max_index = ( sizeof(charset) - 1 );
        return charset[ rand() % max_index ];
    };
    std::string str( UIDLEN, 0 );
    std::generate_n( str.begin(), UIDLEN, randchar );

    return str;
};

treenode::treenode() : uid{0}, children{0}, parent{nullptr}
{
    uid = treenode::generate_random_uid();
    parent = nullptr;
};

treenode::treenode( const std::string& uid, treenode *parent, std::vector<treenode *>& children )
{
    this->uid = uid;
    this->parent = parent;
    this->children = children;
}

treenode::~treenode()
{
    parent = nullptr;
};

void treenode::insert_child( treenode *child )
{
    child->set_parent( this );
    children.insert( children.begin(), child );
};

void treenode::add_child_at( treenode *child, const int& position )
{
    tree_const_it it_pos = children_at( position );
    child->set_parent( this );
    children.insert( it_pos, child );
};

const int& treenode::append_child( treenode *child )
{
    child->set_parent( this );
    children.push_back( child );
    return children.size();
};

const int& treenode::erase_child( treenode *child )
{
    tree_const_it it, ie;
    for ( it = children_begin(); it != ie; it++ ) {
        if ( *it == child ) {
            children.erase( it );
        }
    }
    return it - children_begin();
};

void treenode::pop_child_at( const uint& position )
{
    tree_const_it it = children_at( position ); 
    children.erase( it );
}

treenode* treenode::get_parent()
{
    return parent;
};

void treenode::set_parent( treenode *parent )
{
    this->parent = parent;
}

treenode::tree_const_it treenode::children_at( const uint& position ) {
    return children.begin() + position;
}
