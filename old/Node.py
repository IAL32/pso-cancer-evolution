import random
import networkx as nx
from ete3 import Tree
from graphviz import Source
import operator

import json

class Node(Tree):

    def _get_uid(self):
        return hex(id(self))

    def __init__(self, name, parent, mutation_id, loss=False):
        self.mutation_id = mutation_id
        self.loss = loss

        super().__init__(newick=None, name=name)

        if parent: # automatically add this node to its parent on creation
            parent.add_child(self)

    uid = property(fget=_get_uid)

    def __str__(self):
        return str(self.name) + ("-" if self.loss else "")

    def __repr__(self):
        return str(self.name) + ("-" if self.loss else "")

    def fix_for_losses(self, helper, tree, delete_only=False):
        # saving current children list, it will change if we delete
        # the current node
        if helper.k == 0:
            return
        children = [c for c in self.children]
        for n in children:
            n.fix_for_losses(helper, tree)

        if self.loss and self in tree.losses_list:
            valid = self.is_loss_valid()
            lost = self.is_mutation_already_lost(self.mutation_id, k=helper.k)

            if (not valid) or lost:
                if not delete_only:
                    self.delete_b(helper, tree)
                else:
                    self.delete(prevent_nondicotomic=False)
    
    def delete_b(self, helper, tree):
        tree.losses_list.remove(self)
        tree.k_losses_list[self.mutation_id] -= 1
        self.delete(prevent_nondicotomic=False)

        # TODO: workout how can sigma be done
        # for i in range(helper.cells):

    def find_node_by_uid(self, uid_):
        return next(self.iter_search_nodes(uid=uid_), None)

    def is_loss_valid(self, mutation_id=None):
        """ Checks if current node mutation is valid up until the root node """
        if mutation_id is None:
            mutation_id = self.mutation_id
        for par in self.iter_ancestors():
            if par.mutation_id == mutation_id:
                return True
        return False

    def is_mutation_already_lost(self, mutation_id, k=3):
        """
            Checks if mutation is already lost in the current tree
        """
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

    def get_depth(self):
        ancestors = [n for n in self.iter_ancestors()]
        return len(ancestors)

    def get_height(self):
        "Returns the tree height from the current node"
        height = 0
        for child in self.children:
            height = max(height, child.get_height())
        return height + 1

    def copy_from(self, node):
        self.name = node.name
        self.mutation_id = node.mutation_id
        self.loss = node.loss

    def swap(self, node):
        """ Switch this data with with that of another node """
        tmp_node = Node(self.name, None, self.mutation_id, self.loss)
        self.copy_from(node)
        node.copy_from(tmp_node)

    def get_clade_distance(self, helper, nclades, tree_mt, distance, root=False):
        clade = None
        clade_mut_number = None
        nodes = self.get_clades() if not root else self.get_cached_content()

        for cl in nodes:
            mutations, mut_number = cl.mutation_number(helper)
            if mut_number <= distance and mut_number >= tree_mt - distance:
                if clade is None or mut_number > clade_mut_number:
                    clade = cl
                    clade_mut_number = mut_number
        return clade

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

            This method returns a list of tuples with two elements for each item:
            mutations(T) = { (m, n) : m = mutation_number(n), n € T }
        """

        # sommo il numero di mutazioni acquisite per ogni nodo dell'albero
        nodes = self.get_cached_content()
        mutations = []
        s = 0
        for n in nodes:
            n_genotype = [0] * helper.mutations
            n.get_genotype_profile(n_genotype)
            sum_ = 0
            for m in n_genotype:
                sum_ += m

            mutations.append((sum_, n))
            s += sum_
        return mutations, s

    def get_clades_max_nodes(self, max=1):
        clades = []
        for cl in self.get_clades():
            if len(cl.get_cached_content()) <= max and not cl.loss:
                clades.append(cl)
        return clades

    def distance(self, helper, tree):
        """
            Calculates the distance between this tree and another.
            Tree: if compared with the same tree, it has to be a copy of it,
            not using the same reference, otherwise we will get errors.
            The formula we use in order to calculate the distance between two trees
            is as follows:

            d(T1, T2) = max ( sum_{x € T1}(m(x)), sum_{x € T2}(m(x)) ) - max_weight_matching(x)
            d(T1, T2) € [ 0; (m' * (m' + 1)) / 2 - (m * (m + 1)) / 2 ]

        """
        clades_t1 = self.get_clades()
        clades_t2 = tree.get_clades()

        G = nx.Graph()
        G.add_nodes_from(clades_t1, bipartite=0)
        G.add_nodes_from(clades_t2, bipartite=1)
        edges = []
        weights = []

        mutations_t1, mut_number_t1 = self.mutation_number(helper)
        mutations_t2, mut_number_t2 = tree.mutation_number(helper)

        for cl1 in clades_t1:
            for cl2 in clades_t2:
                # w(e) = n. common mutations between the two clades
                w = Node.common_clades_mutation(helper, cl1, cl2)
                edges.append((cl1, cl2))
                weights.append(w)
                G.add_edge(cl1, cl2, weight=w)
        # max weight matching

        max_matching = nx.algorithms.matching.max_weight_matching(G)

        max_weight = 0
        max_weight_edge = None
        matched_weights = []

        for (kk, vv) in max_matching:
            # searching for that edge position
            for i, (k, v) in enumerate(edges):
                if kk == k and vv == v or kk == v and vv == k:
                    # max_weight += weights[i]
                    matched_weights.append(weights[i])
                    break
            if max_weight == 0 or max_weight < weights[i]:
                max_weight = weights[i]
                max_weight_edge = edges[i]
        distance = max(mut_number_t1, mut_number_t2) - max_weight

        print(matched_weights)
        print(max_matching)
        return distance
        # return distance, max_weight_edge
        # return distance, mutations_t1, mut_number_t1, mutations_t2, mut_number_t2

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
    
    def check_integrity(self):
        for n in self.traverse():
            for c in n.children:
                assert(c.up == n)

    def attach_clade(self, helper, tree, clade):
        "Remove every node already in clade"
        root = self.get_tree_root()
        nodes_list = root.get_cached_content()
        # clade.fix_for_losses(helper, tree, delete_only=True)
        clade_nodes = clade.get_cached_content()
        clade_to_attach = self
        for cln in clade_nodes:
            removed = []
            if cln.loss:
                if clade_to_attach.is_mutation_already_lost(cln.mutation_id):
                    cln.delete(prevent_nondicotomic=False)
                else:
                    tree.losses_list.append(cln)
                    tree.k_losses_list[cln.mutation_id] += 1
            else:
                for n in nodes_list:
                    if n.mutation_id != -1 and cln.mutation_id == n.mutation_id and not n.loss:
                        # moving up
                        if clade_to_attach == n:
                            clade_to_attach = n.up
                        
                        n.delete(prevent_nondicotomic=False)
                        removed.append(n)
            for r in removed:
                nodes_list.pop(r)
        clade_to_attach.add_child(clade)

    def attach_clade_and_fix(self, helper, tree, clade):
        """
        Attaches a clade to the phylogeny tree and fixes everything
        """
        for n in clade.traverse():
            if not tree.k_losses_list[n.mutation_id] <= helper.k:
                n.delete(prevent_nondicotomic=False)
        self.attach_clade(helper, tree, clade)
        self.fix_for_losses(helper, tree)
        # self.check_integrity()

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
            props = {"label": "%s" % (n.name)}
            if n.loss: # marking back-mutations
                props["color"] = "red"
            #     for p in n.iter_ancestors():
            #         if n.mutation_id == p.mutation_id and not p.loss:
            #             out += self._to_dot_node(n.uid, p.uid, props={"style": "dashed", "color": "gray"})
            #             break
            out += self._to_dot_node(n.uid, props=props)
            out += self._to_dot_node(self.uid, n.uid)
            if not n.is_leaf():
                out += n.to_dot()

        if not self.up: # first
            out += '\n}\n'
        return out

    def _to_json_children(self):
        """
            Support function for printing the json tree
        """
        node = {"name": self.name, "uid": self.uid, "loss": self.loss, "children": []}
        for n in self.children:
            node["children"].append(n._to_json_children())
        return node

    def to_json(self):
        """
            Returns a json string representing the current tree
        """
        node = {"name": self.name, "uid": self.uid, "loss": self.loss, "root": True, "children": []}
        for n in self.children:
            node["children"].append(n._to_json_children())
        return json.dumps(node, indent=4)

    def to_string(self):
        return "[uid: %s; dist: %d]" % (str(self.uid), self.get_distance(self.get_tree_root()))

    def _to_tikz_node(self):
        out = ''
        back_mutation = ''
        for c in self.get_children():
            out += c._to_tikz_node()
        if self.loss:
            back_mutation = ',color=red'
        return '\n\t[{%s}%s %s]' % (self.name, back_mutation, out)

    def to_tikz(self):
        nodes = self.get_cached_content()
        # refer to official "forest" package documentation for forked edges
        # changed it a bit:
        # \forestset{
        # 	declare dimen={fork sep}{0.5em},
        # 	forked edge'/.style={
        # 		edge={rotate/.option=!parent.grow},
        # 		edge path'={(!u.parent anchor) -- ++(\forestoption{fork sep},0) |- (.child anchor)},
        # 	},
        # 	forked edge/.style={
        # 		on invalid={fake}{!parent.parent anchor=children},
        # 		child anchor=parent,
        # 		forked edge',
        # 	},
        # 	forked edges/.style={for nodewalk={#1}{forked edge}},
        # 	forked edges/.default=tree,
        # 	aligned terminal/.style={if n children=0{
        # 		tier=terminal
        # 	}{}},
        # 	germline/.style={
        # 		for tree = {grow'=0,draw,aligned terminal}, forked edges
        # 	}
        # }
        out = '\\begin{forest}\n\tgermline'
        out += '\n\t[{%s} ' % self.name
        for c in self.get_children():
            out += c._to_tikz_node()
        return out + ']\n\\end{forest}'

    def save(self, filename="test.gv", fileformat="dot"):
        if fileformat == "dot":
            Source(self.to_dot(), filename=filename, format="png").render()
        elif fileformat == "json":
            with open(filename, 'w') as f:
                f.write(self.to_json())