
class Data(object):

    def __init__(self, nofparticles, iterations):
        self.pso_start = 0
        self.pso_end = 0
        self.initialization_start = 0
        self.initialization_end = 0
        self.initialization_times = []
        self.starting_likelihood = 0

        self.nofparticles = nofparticles
        self.iterations = iterations
        self.particle_iteration_times = [[] for p in range(nofparticles)]
        self.iteration_times = []

        self.iteration_new_particle_best = [[0 for p in range(nofparticles)] for n in range(iterations)]
        self.iteration_new_best = [[0 for p in range(nofparticles)] for n in range(iterations)]

        self.false_negatives = 0
        self.false_positives = 0
        self.missing_values = 0

    def pso_passed_seconds(self):
        return self._passed_seconds(self.pso_start, self.pso_end)

    def initialization_passed_seconds(self):
        return self._passed_seconds(self.initialization_start, self.initialization_end)

    def iteration_passed_seconds(self, start, end):
        return self._passed_seconds(start, end)

    def average_iteration_time(self):
        """
            Returns the average time for every iteration
        """
        sum = 0
        for t in self.iteration_times:
            sum += t
        return t / self.iterations

    def average_particle_time(self):
        """
            For each particle, return the average time
            that particle has taken to complete its step
        """
        sums = [0] * self.nofparticles

        for i in range(self.nofparticles):
            for j in range(self.iterations):
                sums[i] += self.particle_iteration_times[i][j]

        averages = [0] * self.nofparticles
        for p, s in enumerate(sums):
            averages[p] = s / self.iterations
        
        return averages

    def average_overall_particle(self):
        p_times = self.average_particle_time()

        sum = 0
        for t in p_times:
            sum += t
        
        return sum / self.nofparticles

    def _passed_seconds(self, start, end):
        return end - start
