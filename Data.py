from Tree import Tree
import matplotlib.pyplot as plt
import random
import sys
class Data(object):

    def __init__(self, nofparticles, iterations, seed):
        self.pso_start = 0
        self.pso_end = 0
        self.initialization_start = 0
        self.initialization_end = 0
        self.initialization_times = []
        self.starting_likelihood = 0
        if seed == -1:
            seed = random.randrange(sys.maxsize)
        self.seed = seed

        random.seed(seed)

        self.nofparticles = nofparticles
        self.iterations = iterations
        self.particle_iteration_times = [[] for p in range(nofparticles)]
        self.iteration_times = []

        self.iteration_new_particle_best = [[0 for p in range(nofparticles)] for n in range(iterations)]
        self.iteration_new_best = [[0 for p in range(nofparticles)] for n in range(iterations)]

        self.true_positive = 0
        self.true_negative = 0
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

    def average_iteration_particle_time(self):
        average_times = [0] * self.iterations
        for j in range(self.iterations):
            for i in range(self.nofparticles):
                average_times[j] += self.particle_iteration_times[i][j]
        
        for j in range(self.iterations):
            average_times[j] /= self.iterations
        
        return average_times

    def average_overall_particle(self):
        p_times = self.average_particle_time()

        sum = 0
        for t in p_times:
            sum += t

        return sum / self.nofparticles

    def _passed_seconds(self, start, end):
        return end - start

    def summary(self, helper):
        Tree.greedy_loglikelihood(helper, helper.best_particle.best, self)
        print ("Number of particles: %d" % self.nofparticles)
        print ("Seed used: %d" % self.seed)
        print ("Number of iterations: %d" % self.iterations)
        print ("Number of cells: %d" % helper.cells)
        print ("Number of mutations: %d" % helper.mutations)
        print ("Starting likelihood: %f" % self.starting_likelihood)
        print ("Best likelihood: %f" % helper.best_particle.best.likelihood)
        print ("Added mutations:", helper.best_particle.best.losses_list)
        print ("False negatives: %d" % self.false_negatives)
        print ("False positives: %d" % self.false_positives)
        print ("Added missing values: %d" % self.missing_values)
        print ("PSO completed in %f seconds" % (self.pso_passed_seconds()))
        print ("Initialization took %f seconds" % self.initialization_passed_seconds())
        print ("Average iteration time: %f seconds" % self.average_iteration_time())
        print ("Average particle time: %f seconds" % self.average_overall_particle())

        average_particle_time_ = self.average_particle_time()
        plt.plot(average_particle_time_)
        plt.ylim(bottom=0, top=max(average_particle_time_))
        plt.savefig("results/average_particle_time.png")

        plt.clf()
        average_iteration_particle_time_ = self.average_iteration_particle_time()
        plt.plot(average_iteration_particle_time_)
        plt.ylim(top = max(average_iteration_particle_time_))
        plt.savefig("results/average_iteration_particle_time.png")

        helper.best_particle.best.phylogeny.save("results/best.gv")