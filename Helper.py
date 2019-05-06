class Helper(object):
    def __init__(self, matrix, mutations, mutation_names, cells, alpha, beta, k):
        self.matrix = matrix
        self.mutations = mutations
        self.mutation_names = mutation_names
        self.cells = cells
        self.alpha = alpha
        self.beta = beta
        self.k = k
        self.best_particle = None
