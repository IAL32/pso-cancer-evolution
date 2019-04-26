class Particle(object):

    def __init__(self, tree, best_likelihood, best_sigma, losses_list, k_losses_list, velocity, operations, original):
        self.swap(tree, best_likelihood, best_sigma, losses_list, k_losses_list, velocity, operations, original)

    def swap(self, tree, best_likelihood, best_sigma, losses_list, k_losses_list, velocity, operations = [], original = None):
        self.tree = tree
        self.best_likelihood = best_likelihood
        self.best_sigma = best_sigma
        self.losses_list = losses_list
        self.k_losses_list = k_losses_list
        self.velocity = velocity

        self.operations = operations
        self.original = original

    @classmethod
    def fromscratch(cls, cells, mutations, tree):
        best_likelihood = float("-inf")
        best_sigma = [None] * cells
        losses_list = []
        k_losses_list = [0] * mutations
        velocity = 0
        operations = []
        original = tree.copy()
        return cls(tree, best_likelihood, best_sigma, losses_list, k_losses_list, velocity, operations, original)
    
    def __repr__(self):
        return "bl: " + str(self.best_likelihood)