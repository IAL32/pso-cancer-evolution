class Helper(object):
    def __init__(self, matrix, mutations, mutation_names, cells, alpha, beta, k, c1, c2, inertia, vmax):
        self.matrix = matrix
        self.mutations = mutations
        self.mutation_names = mutation_names
        self.cells = cells
        self.alpha = alpha
        self.beta = beta
        self.k = k
        self.c1 = c1
        self.c2 = c2
        self.inertia = inertia
        self.vmax = vmax
        self.best_likelihood = float("-inf")
        self.best_likelihood_particle = None