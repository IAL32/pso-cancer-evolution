import copy
from Tree import Tree

class Particle(object):

    def _last_tree_get(self):
        return self.trees[-1]

    def __init__(self, cells, mutations, mutation_names):
        # tree linked list
        self.trees = [Tree.random(cells, mutations, mutation_names)]
        # best tree with the best likelihood so far
        self.best = self.trees[-1]

    last_tree = property(fget=_last_tree_get)

    def __repr__(self):
        return "bl: " + str(self.best.likelihood)