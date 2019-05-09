import random
import string

import matplotlib.pyplot as plt
import networkx as nx
from ete3 import Tree
from graphviz import Source

random.seed(1)

def rid(k=6):
    return ''.join(random.choices(string.digits + 'abcdef', k=k))
class Node(Tree):

    def __init__(self, name, parent, mutation_id, uid, loss=False):
        self._up = parent
        self.uid = uid
        self.mutation_id = mutation_id
        self.loss = loss

        super().__init__(newick=None, name=name)

        if parent: # automatically add this node to its parent on creation
            parent.add_child(self)

    def __str__(self):
        return self.name + ("-" if self.loss else "")

    def fix_for_losses(self, helper, tree):
        # saving current children list, it will change if we delete
        # the current node
        children = [c for c in self.children]
        for n in children:
            n.fix_for_losses(helper, tree)

        if self.loss and self in tree.losses_list:
            valid = self.is_loss_valid()
            lost = self.is_mutation_already_lost(self.mutation_id)

            if (not valid) or lost:
                self.delete_b(helper, tree)
    
    def delete_b(self, helper, tree):
        tree.losses_list.remove(self)
        tree.k_losses_list[self.mutation_id] -= 1
        self.delete(prevent_nondicotomic=False)

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

    def swap(self, node):
        """ Switch this data with with that of another node """
        tmp_node = Node(self.name, None, self.mutation_id, self.uid, self.loss)
        self.copy_from(node)
        node.copy_from(tmp_node)

    def get_clades_at_height(self, height=1):
        " Returns a list of clades at the desired height "
        clades = self.get_clades()
        for cl in clades:
            if cl.get_height() != height:
                clades.remove(cl)
        return clades

    def _get_parent_at_height(self, height=1):
        " Support function that returns the parent node at the desired height "

        par = self.up
        climb = 0
        while (par is not None and climb < height):
            climb += 1
            par = par.up
        return par

    def get_clades(self):
        """
        Clades are defined as every node in the tree, excluding the root
        """
        if self.mutation_id != -1:
            raise SystemError("Cannot get clades from a non-root node!")
        nodes_list = list(self.get_cached_content().keys())
        nodes_list.remove(self)
        return nodes_list

    def get_genotype_profile(self, genotypes):
        " Walks up to the root and maps the genotype for the current node mutation "
        if self.mutation_id == -1:
            return
        if not self.loss:
            genotypes[self.mutation_id] += 1
        else:
            genotypes[self.mutation_id] -= 1

        self.up.get_genotype_profile(genotypes)

    def get_random_node(self):
        "Returns a random node"
        return random.choice(list(self.get_cached_content().keys()))

    def mutation_number(self, helper):
        """
            Suppose that we have the following tree:
            T:
                       /-c
                    /d|
            -germline  \-a
                   |
                    \-b
            
            Our tree can be represented by the following matrix,
            obtained by combining every mutation genotype:
            M(T) =
            a = 1 0 0 1
            b = 0 1 0 0
            c = 0 0 1 1
            d = 0 0 0 1
            And the sum of mutations is the sum of very 1 in the matrix:
            a = 1 + 0 + 0 + 1 = 2
            b = 0 + 1 + 0 + 0 = 1
            c = 0 + 0 + 1 + 1 = 2
            d = 0 + 0 + 0 + 1 = 1
            a + b + c + d     = 6
            And 6 is the sum of the number of mutations in the tree.
        """
        sum = 0
        # sommo il numero di mutazioni acquisite per ogni nodo dell'albero
        for n in self.get_cached_content():
            n_genotype = [0] * helper.mutations
            n.get_genotype_profile(n_genotype)
            for m in n_genotype:
                sum += m
        return sum

    def distance(self, helper, tree):
        """
            Calculates the distance between this tree and another.
            Tree: if compared with the same tree, it has to be a copy of it,
            not using the same reference, otherwise we will get errors.
            The formula we use in order to calculate the distance between two trees
            is as follows:

            d(T1, T2) = max ( sum_{x € T1}(m(x)), sum_{x € T2}(m(x)) ) - max_weight_matching(x)
            d(T1, T2) € [ m; (m * (m + 1)) / 2 ]

        """
        clades_t1 = self.get_clades()
        clades_t2 = tree.get_clades()

        G = nx.Graph()
        G.add_nodes_from(clades_t1, bipartite=0)
        G.add_nodes_from(clades_t2, bipartite=1)
        edges = []
        weights = []

        mutations_t1 = self.mutation_number(helper)
        mutations_t2 = tree.mutation_number(helper)

        for cl1 in clades_t1:
            for cl2 in clades_t2:
                # w(e) = n. common mutations between the two clades
                w = Node.common_clades_mutation(helper, cl1, cl2)
                edges.append((cl1, cl2))
                weights.append(w)
                G.add_edge(cl1, cl2, weight=w)
        # max weight matching

        max_matching = list(nx.algorithms.matching.max_weight_matching(G))
        max_weight = 0
        for (kk, vv) in max_matching:
            # searching for that edge position
            for i, (k, v) in enumerate(edges):
                if kk == k and vv == v or kk == v and vv == k:
                    max_weight += weights[i]
                    break
        print (max_matching)
        distance = max(mutations_t1, mutations_t2) - max_weight

        return distance

    def back_mutation_ancestry(self):
        """
            Returns a list of nodes representing where a back mutation
            happened. Mostly used to know where NOT to cut.
        """
        back_mutations = []
        for p in self.iter_ancestors():
            if p.loss:
                back_mutations.append(p)
        return back_mutations

    def attach_clade(self, clade):
        "Remove every node already in clade"
        nodes_list = list(self.get_cached_content().keys())
        clade_nodes_list = clade.get_cached_content()
        for cln in clade_nodes_list:
            for n in nodes_list:
                if n.mutation_id != -1 and cln.name == n.name:
                    n.delete(prevent_nondicotomic=False)
                    nodes_list.remove(n)

        self.add_child(clade)

    def attach_clade_and_fix(self, helper, tree, clade):
        """
        Attaches a clade to the phylogeny tree and fixes everything
        """
        self.attach_clade(clade)

        losses_list, k_losses_list = tree.calculate_losses_list(helper.k)
        tree.losses_list = losses_list
        tree.k_losses_list = k_losses_list
        self.fix_for_losses(helper, tree)

    @classmethod
    def common_clades_mutation(cls, helper, clade1, clade2):
        """
            This function can be seen as the logic and between
            two binary strings, and then the sum between every element.
            Suppose we have the following trees:
            T1:
                       /-c
                    /d|
            -germline  \-a
                   |
                    \-b
            T2:
                       /-d
                    /c|
            -germline  \-b
                   |
                    \-a
            
            And suppose we are comparing the clades c1 and b2.
            genotype(c1) = [0 0 1 1]
            genotype(b2) = [0 1 1 0]
            logic_and = [0 0 1 1] & [0 1 1 0] = [0 0 1 0]
            sum = 0 + 0 + 1 + 0 = 1
        """

        clade1_genotype = [0] * helper.mutations
        clade2_genotype = [0] * helper.mutations
        common = 0

        # ignoring back mutations
        clade1.get_genotype_profile(clade1_genotype)
        clade2.get_genotype_profile(clade2_genotype)

        for m in range(helper.mutations):
            if clade1_genotype[m] == clade2_genotype[m] == 1:
                common += 1
        return common

    def _to_dot_label(self, d={}):
        """

        Returns a string representing the list of properties
        indicated by d.
        Ex.: d = {
            "label": "name",
            "color": "red"
        }
        Will result in:
        [label="name",color="red"]
        """
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

    def to_dot(self, root=False):
        out = ''
        if not self.up or root: # first graph node
            out += 'graph {\n\trankdir=UD;\n\tsplines=line;\n\tnode [shape=circle]'
            out += self._to_dot_node(self.uid, props={"label": self.name})
        for n in self.children:
            props = {"label": "%s\nuid: %s" % (n.name, str(n.uid))}
            if n.loss: # marking back-mutations
                props["color"] = "red"
                for p in n.iter_ancestors():
                    if n.mutation_id == p.mutation_id and not p.loss:
                        out += self._to_dot_node(n.uid, p.uid, props={"style": "dashed", "color": "gray"})
                        break
            out += self._to_dot_node(n.uid, props=props)
            out += self._to_dot_node(self.uid, n.uid)
            if not n.is_leaf():
                out += n.to_dot()

        if not self.up: # first
            out += '\n}\n'
        return out

    def to_string(self):
        return "[uid: %s; dist: %d]" % (str(self.uid), self.get_distance(self.get_tree_root()))

    def save(self, filename="test.gv"):
        Source(self.to_dot(), filename=filename, format="png").render()
