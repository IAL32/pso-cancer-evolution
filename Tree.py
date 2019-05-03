from Node import Node, rid
import random as r
import copy
r.seed(1)

class Tree(object):

    def __init__(self, cells, mutations):
        self.cells = cells
        self.mutations = mutations
        self.losses_list = []
        self.k_losses_list = [0] * mutations
        self.best_sigma = [0] * cells
        self.likelihood = float("-inf")
        self.phylogeny = None

    def calculate_losses_list(self, mutations, k):
        losses_list = []
        k_losses_list = [0] * mutations
        for n in self.phylogeny.traverse():
            if n.loss:
                losses_list.append(n)
                k_losses_list[n.mutation_id] += 1

                # sanity check
                if k_losses_list[n.mutation_id] > k:
                    raise SystemError("Mutation %d has too many losses!" % n.mutation_id)
        return losses_list, k_losses_list

    def copy(self):
        "Copies everything in this tree"
        t = Tree(self.cells, self.mutations)
        t.likelihood = self.likelihood
        t.phylogeny = self.phylogeny.copy()
        for n in t.phylogeny.traverse():
            if n.loss:
                t.losses_list.append(n)
        t.k_losses_list = [l for l in self.k_losses_list]
        return t

    @classmethod
    def random(cls, cells, mutations, mutation_names):
        t = Tree(cells, mutations)
        random_tree = cls._generate_random_btree(mutations, mutation_names)
        t.phylogeny = random_tree
        return t

    @classmethod
    def _generate_random_btree(cls, mutations, mutation_names):
        """ Generates a random binary tree """
        root = Node("germline", None, -1, 0)

        rantree = [i for i in range(mutations)]

        r.shuffle(rantree, r.random)

        nodes = [root]
        append_node = 0
        i = 0
        while i < mutations:
            nodes.append(
                Node(mutation_names[rantree[i]], nodes[append_node], rantree[i], rid())
            )
            i += 1

            if i < mutations:
                nodes.append(
                    Node(mutation_names[rantree[i]], nodes[append_node], rantree[i], rid())
                )
            append_node += 1
            i += 1
        
        return root