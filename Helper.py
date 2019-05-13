class Helper(object):
    def __init__(self, matrix, mutations, mutation_names, cells, alpha, beta, k, c1, c2):
        self.matrix = matrix
        self.mutations = mutations
        self.mutation_names = mutation_names
        self.cells = cells
        self.alpha = alpha
        self.beta = beta
        self.k = k
        self.c1 = c1
        self.c2 = c2
        self.best_particle = None

    def max_phylogeny_mutations(self):
        """
            As of now, we are using the sum of the number of mutations of the
            phylogeny in order to calculate the distance.
            This number is can be between the range: [m; (m(m + 1)) / 2]
            The first is given in the case we have the following case:
                /-a
               |
               |--b
            -root
               |--c
               |
                \-d
            Where every node has at most one mutation, hence the matrix representing this tree
            will be:
            a = 1 0 0 0
            b = 0 1 0 0
            c = 0 0 1 0
            d = 0 0 0 1
            Which has the sum of the mutations = m
            The extreme case is where we have all the mutations aligned in a single branch:
            -root -- a -- b -- c -- d
            a = 1 0 0 0
            b = 1 1 0 0
            c = 1 1 1 0
            d = 1 1 1 1
            Which has the sum of the mutations = (m(m + 1)) / 2
        """
        return (self.mutations * (self.mutations + 1)) / 2
