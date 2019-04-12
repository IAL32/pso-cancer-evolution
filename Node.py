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

    def fix_for_losses(self, helper, tree_helper):

        if not self.is_leaf():
            self.children[0].fix_for_losses(helper, tree_helper)
        next_sibling = self.next_sibling()
        if next_sibling:
            next_sibling.fix_for_losses(helper, tree_helper)

        if self.loss:
            valid = self.is_loss_valid()
            lost = self.is_mutation_already_lost(self.mutation_id)

            if not valid or lost:
                self.delete_b(helper, tree_helper)
    
    def delete_b(self, helper, tree_helper):
        tree_helper.losses_list.remove(self)
        tree_helper.k_losses_list[self.mutation_id] -= 1
        self.delete()

        # TODO: workout how can sigma be done
        # for i in range(helper.cells):

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
        return (self in node.iter_ancestors())

    def prune_and_reattach(self, node_reattach):
        """ Detaches current node (with all its descendants) and reattaches it into another node """
        self.detach()
        node_reattach.add_child(self)
    
    def copy_from(self, node):
        self.uid = node.uid
        self.name = node.name
        self.mutation_id = node.mutation_id
        self.loss = node.loss

    def replace(self, node):
        """ Switch this data with with that of another node """
        tmp_uid = self.uid
        tmp_name = self.name
        tmp_mutation_id  = self.mutation_id
        tmp_loss = self.loss

        self.copy_from(node)

        node.uid = tmp_uid
        node.name = tmp_name
        node.mutation_id = tmp_mutation_id
        node.loss = tmp_loss
    
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
            props = {"label": n.name + "\nL: " + str(self.loss)}
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

    def save(self, filename="test.gv"):
        Source(self.to_dot(), filename=filename, format="png").render()