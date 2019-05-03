from ete3 import Tree
import random
import string
from graphviz import Source

def rid(k=6):
    return ''.join(random.choices(string.digits + 'abcdef', k=k))
class Node(Tree):

    def _get_uid(self):
        return self._uid
    def _set_uid(self, uid):
        self._uid = uid

    def _get_mutation_id(self):
        return self._mut_id
    def _set_mutation_id(self, mutation_id):
        self._mut_id = mutation_id

    def _get_loss(self):
        return self._loss
    def _set_loss(self, loss):
        self._loss = loss

    def __init__(self, name, parent, mutation_id, uid, loss=False):
        self._up = parent
        self._uid = uid
        self._mut_id = mutation_id
        self._loss = loss

        super().__init__(newick=None, name=name)

        if parent: # automatically add this node to its parent on creation
            parent.add_child(self)

    uid = property(fget=_get_uid, fset=_set_uid)
    mutation_id = property(fget=_get_mutation_id, fset=_set_mutation_id)
    loss = property(fget=_get_loss, fset=_set_loss)

    def fix_for_losses(self, helper, tree):

        for n in self.traverse():
            if n.uid != self.uid:
                n.fix_for_losses(helper, tree)

        if self.loss:
            valid = self.is_loss_valid()
            lost = self.is_mutation_already_lost(self.mutation_id)

            if (not valid) or lost:
                self.delete_b(helper, tree)
    
    def delete_b(self, helper, tree):
        tree.losses_list.remove(self)
        tree.k_losses_list[self.mutation_id] -= 1
        self.delete()

        # TODO: workout how can sigma be done
        # for i in range(helper.cells):

    def find_node_by_uid(self, uid_):
        return next(self.iter_search_nodes(uid=uid_), None)

    def copy_all(self):
        cached_content = self.get_cached_content(leaves_only=False)
        return self.copy_until_depth(len(cached_content))
    def copy_until_depth(self, depth):
        # root node
        new_node = Node(self.name, None, self.mutation_id, self.uid, self.loss)
        self._copy_until_depth(depth - 1, new_node)
        return new_node
    def _copy_until_depth(self, depth, up):
        if depth == 0:
            return
        for c in self.children:
            new_child = Node(c.name, None, c.mutation_id, c.uid, c.loss)
            up.add_child(new_child)
            c._copy_until_depth(depth - 1, new_child)

    def is_loss_valid(self):
        """ Checks if current node mutation is valid up until the root node """
        for par in self.iter_ancestors():
            if par.mutation_id == self.mutation_id:
                return True
        return False

    def is_mutation_already_lost(self, mutation_id):
        """ Checks if mutation is already lost in the current tree """
        for par in self.iter_ancestors():
            if par.loss and par.mutation_id == mutation_id:
                return True
        
        return False

    def previous_sibling(self):
        " Returns the previous sibling, if it exists "
        if self.up is None:
            return None
        if len(self.up.children) == 1:
            return None
        self_position = self.up.children.index(self)
        if self_position - 1 == -1:
            return None
        sibling_position = self_position - 1
        return self.up.children[sibling_position]

    def next_sibling(self):
        " Returns the next sibling, if it exists "
        if self.up is None:
            return None
        if len(self.up.children) == 1: # single node
            return None
        self_position = self.up.children.index(self)
        if self_position + 1 == len(self.up.children):
            return None
        sibling_position = self_position + 1
        return self.up.children[sibling_position]

    def is_ancestor_of(self, node):
        """ Checks if current node is parent of the given arguent node """
        par = self.up
        while par != None:
            if par.uid == node.uid:
                return True
            par = par.up
        return False

    def prune_and_reattach(self, node_reattach):
        """ Detaches current node (with all its descendants) and reattaches it into another node """

        if node_reattach.is_ancestor_of(self):
            return 1
        if node_reattach.up.uid == self.uid:
            return 1
        if self.up is None:
            return 1
        if self.uid == node_reattach.uid:
            return 1

        self.detach()
        node_reattach.add_child(self)
        return 0
    
    def get_height(self):
        "Returns the tree height from the current node"
        
        height = 0
        for child in self.children:
            height = max(height, child.get_height())
        return height + 1
    def copy_from(self, node):
        self.uid = node.uid
        self.name = node.name
        self.mutation_id = node.mutation_id
        self.loss = node.loss

    def swap(self, node, tree):
        """ Switch this data with with that of another node """
        if self.loss:
            tree.losses_list.remove(self)
        tmp_uid = self.uid
        tmp_name = self.name
        tmp_mutation_id  = self.mutation_id
        tmp_loss = self.loss

        self.copy_from(node)

        node.uid = tmp_uid
        node.name = tmp_name
        node.mutation_id = tmp_mutation_id
        node.loss = tmp_loss

        if node.loss:
            tree.losses_list.append(node)

    def get_clades_at_height(self, height=1):
        " Returns a list of clades at the desired height "
        tmp = []
        # getting only leaves from cache, faster than iterating
        # over them
        cached_nodes = self.get_cached_content()

        for n in cached_nodes:
            if n.is_leaf():
                for cl in n._get_parent_at_height(height):
                    if len(tmp) > 0:
                        if cl not in tmp:
                            # temporary adding clade, as it will be checked later anyways
                            tmp.append(cl)
                    else:
                        tmp.append(cl)
        clades = []

        # check that every candidate clade is not an ancestor of someone else,
        # otherwise we could get an overlap of clades

        for cl in tmp:
            add = True
            for cl_ in tmp:
                # skip same node
                if not cl.uid == cl_.uid:
                    # check all possible clades
                    if cl_.is_ancestor_of(cl):
                        add = False
                        break
            if add:
                clades.append(cl)

        return clades

    def _get_parent_at_height(self, height=1):
        " Support function that yields a clade at the desired height "
        if not self.is_leaf():
            raise SystemError("Trying to get clades from a non-leaf node")

        par = self.up
        climb = 0
        while (par is not None and climb < height):
            yield par
            climb += 1
            par = par.up

    def _to_dot_label(self, d={}):
        if not len(d):
            return ''

        out = '['
        for i, (key, value) in enumerate(d.items()):
            if isinstance(value, (int, float, complex)):
                out += '%s=%s' % (key, str(value))
            else:
                out += '%s="%s"' % (key, str(value))
            if i < len(d) - 1: # last
                out += ','
        out += ']'
        return out

    def _to_dot_node(self, nodeFromId, nodeToId=None, props={}):
        if nodeToId:
            return '\n\t"%s" -- "%s" %s;' % (nodeFromId, nodeToId, self._to_dot_label(props))
        else: # printing out single node
            return '\n\t"%s" %s;' % (nodeFromId, self._to_dot_label(props))

    def to_dot(self):
        out = ''
        if not self.up: # first graph node
            out += 'graph {\n\trankdir=UD;\n\tsplines=line;\n\tnode [shape=circle]'
            out += self._to_dot_node(self.uid, props={"label": self.name})
        for n in self.children:
            props = {"label": n.name + "\nuid: " + n.uid}
            if n.loss: # marking back-mutations
                props["color"] = "red"
            out += self._to_dot_node(n.uid, props=props)
            out += self._to_dot_node(self.uid, n.uid)
            if not n.is_leaf():
                out += n.to_dot()

        if not self.up: # first
            out += '\n}\n'
        return out

    def get_genotype_profile(self, genotypes):
        " Walks up to the root and maps the genotype for the current node mutation "
        if self.mutation_id == -1:
            return
        if not self.loss:
            genotypes[self.mutation_id] += 1
        else:
            genotypes[self.mutation_id] -= 1

        self.up.get_genotype_profile(genotypes)

    def to_string(self):
        return "[uid: " + str(self.uid) + "; dist: " + str(self.get_distance(self.get_tree_root())) + "]"

    def save(self, filename="test.gv"):
        Source(self.to_dot(), filename=filename, format="png").render()