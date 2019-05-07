import copy
from Operation import Operation
from Tree import Tree

class Particle(object):

    def last_tree(self):
        return self.trees[-1]

    def __init__(self, cells, mutations, mutation_names, number):
        # tree linked list
        self.trees = [Tree.random(cells, mutations, mutation_names)]
        self.number = number
        # best tree with the best likelihood so far
        self.best = self.trees[-1]

    def __repr__(self):
        return "bl: " + str(self.best.likelihood)