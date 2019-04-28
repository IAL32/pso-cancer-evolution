import copy
from Node import Node, rid
import random as r

class Particle(object):

    def __init__(self, tree, best_tree, best_likelihood, best_sigma, losses_list, k_losses_list):
        self.swap(tree, best_tree, best_likelihood, best_sigma, losses_list, k_losses_list)

    def swap(self, tree, best_tree, best_likelihood, best_sigma, losses_list, k_losses_list):
        self.tree = tree
        self.best_tree = best_tree
        self.best_likelihood = best_likelihood
        self.best_sigma = best_sigma
        self.losses_list = losses_list
        self.k_losses_list = k_losses_list

    def calculate_losses_list(self, helper):
        losses_list = []
        k_losses_list = [0] * helper.mutations
        for n in self.tree.traverse():
            if n.loss:
                losses_list.append(n)
                k_losses_list[n.mutation_id] += 1
                
                if k_losses_list[n.mutation_id] > helper.k:
                    raise SystemError("Mutation %d has too many losses!" % n.mutation_id)
        return losses_list, k_losses_list

    def copy(self):
        return Particle(self.tree.copy(), self.best_tree, self.best_likelihood, copy.deepcopy(self.best_sigma), copy.deepcopy(self.losses_list), copy.deepcopy(self.k_losses_list))

    @classmethod
    def fromscratch(cls, cells, mutations, tree):
        best_likelihood = float("-inf")
        best_sigma = [None] * cells
        losses_list = []
        k_losses_list = [0] * mutations

        return cls(tree, tree, best_likelihood, best_sigma, losses_list, k_losses_list)

    @classmethod
    def random(cls, helper):
        return Particle.fromscratch(helper.cells, helper.mutations, Particle._generate_random_btree(helper.mutations, helper.mutation_names))

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


    def __repr__(self):
        return "bl: " + str(self.best_likelihood)