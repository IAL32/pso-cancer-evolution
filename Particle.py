import copy
from Operation import Operation
from Tree import Tree

class Particle(object):

    def __init__(self, cells, mutations, mutation_names, number):
        # tree linked list
        self.current_tree = Tree.random(cells, mutations, mutation_names)
        self.number = number
        # best tree with the best likelihood so far
        self.best = self.current_tree
        self.climb_probability = 1.0

    

    def __repr__(self):
        return "bl: " + str(self.best.likelihood)
