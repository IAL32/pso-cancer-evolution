import multiprocessing as mp
import random
import sys
import time

from Helper import Helper
from Node import Node, rid
from Operation import Operation as Op
from Particle import Particle
from Tree import Tree

random.seed(1)

# global scope for multiprocessing
particles = []
helper = None

def init(nparticles, iterations, matrix, mutations, mutation_names, cells, alpha, beta, k):
    global helper
    global particles
    helper = Helper(matrix, mutations, mutation_names, cells, alpha, beta, k)

    pso(nparticles, iterations, matrix)

    # for i, p in enumerate(particles):
    #     print ("Particle n. %d" % i)
    #     print ("- loglh: %d" % p.best.likelihood)
    #     # p.best.phylogeny.save("trees/tree_last_%d.gv" % i)
    
    print ("Overall best particle: lh: %f" % helper.best_particle.best.likelihood)
    helper.best_particle.best.phylogeny.save("trees/best.gv")

def cb_init_particle(result):
    i, particle = result
    particles[i] = particle
    if (particle.current_tree.likelihood > helper.best_particle.best.likelihood):
        helper.best_particle = particle

def init_particle(i, p, helper):
    lh = Tree.greedy_loglikelihood(helper, p.current_tree)
    p.current_tree.likelihood = lh
    return i, p

def cb_particle_iteration(r):
    result, op, p, tree_copy = r
    # updating log likelihood and bests
    particles[p.number] = p

    if result == 0:
        lh = Tree.greedy_loglikelihood(helper, tree_copy)
        tree_copy.likelihood = lh
        p.current_tree = tree_copy
        # print("Operation %d" % (tree_copy.operation.type))
        if lh != helper.best_particle.best.likelihood and lh > p.best.likelihood:
            # updating particle best
            decreased = (p.best.likelihood - lh) / lh * 100
            print("- !! %d new particle best, before: %d, now: %d, increased by %f%%" % (p.number, p.best.likelihood, lh, decreased))
            p.best = tree_copy
        if lh > helper.best_particle.best.likelihood:
            # updating swarm best
            decreased = (helper.best_particle.best.likelihood - lh) / lh * 100
            print("- !!!!! %d new swarm best, before: %d, now: %d, increased by %f%%" % (p.number, helper.best_particle.best.likelihood, lh, decreased))
            helper.best_particle = p

def particle_iteration(p, helper):
    # print ("Particle n. %d" % i)
    # print ("- loglh: %d" % p.best.likelihood)
    ops = list(range(0, Op.NUMBER))
    result = -1

    while len(ops) > 0 and result != 0:
        op = ops.pop(random.randint(0, len(ops) - 1))

        best_swarm_copy = helper.best_particle.best.copy()
        best_particle_copy = p.best.copy()

        """
        we're going to "mix" three trees:
        - the current tree
        The current solution that has to be modified
        - the best swarm tree
        This guides us to the best solution we know of
        - the best particle tree
        This slows down the way we get to the best solution, in order
        to avoid getting stuck in a local optimum

        We "move" by copying clades from a tree into another at a given height.

        We calculate the distance between:
        current_tree - best_swarm_tree
        current_tree - best_particle_tree

        Given the formula for the distance between phylogenies:
        d(T1, T2) = max ( sum_{x € T1}(m(x)), sum_{x € T2}(m(x)) ) - max_weight_matching(x)

        The more distant we are from an optimum, the faster we have to move.
        The less distant we are from an optimu, the slower we have to move.

        How fast we move defines how high the clades we copy are, so we define it
        as the velocity of our particle. The velocity is directly proportional to
        how distant we are from the other tree.
        v = d(T1, T2)

        Also, we will randomly pick the number "n" of clades, n € [1, 3],
        that we will pick from each tree.
        """

        distance_swarm = p.current_tree.phylogeny.distance(helper, best_swarm_copy.phylogeny)
        distance_particle = p.current_tree.phylogeny.distance(helper, best_particle_copy.phylogeny)

        print(distance_swarm, distance_particle)

        # best_swarm_copy.get_clades_at_height()
        # clades_current_tree = p.current_tree.phylogeny.get_clades()

        # tree_copy = Tree.random(helper.cells, helper.mutations, helper.mutation_names)
        # tree_copy.phylogeny.attach_clade_and_fix(helper, tree_copy, lowest_distance_clade.copy())
        result = Op.tree_operation(helper, tree_copy, op)

    return result, op, p, tree_copy

def pso(nparticles, iterations, matrix):
    global particles
    global helper
    # Particle initialization
    print("Particle initialization...")
    # Random position, each tree is a binary tree at the beginning
    particles = [Particle(helper.cells, helper.mutations, helper.mutation_names, n) for n in range(nparticles)]

    helper.best_particle = particles[0]
    pool = mp.Pool(mp.cpu_count())

    start = time.time()
    processes = []
    # parallelizing tree initialization
    for i, p in enumerate(particles):
        processes.append(pool.apply_async(init_particle, args=(i, p, helper), callback=cb_init_particle))

    for p in processes:
        p.get()
    end = time.time()

    print("Initialization of %d particles: %f seconds" % (len(particles), end - start))

    print ("best logarithm: %f" % helper.best_particle.best.likelihood)

    pool.close()
    pool.join()

    # exit(1)

    # start = time.time()
    # for p in particles:
    #     lh = greedy_tree_loglikelihood(helper, p.current_tree)
    #     p.current_tree.likelihood = lh
    #     if (p.current_tree.likelihood > helper.best_particle.best.likelihood):
    #         helper.best = p
    # end = time.time()
    # print("Standard time for initialization: %f" % (start - end))

    start_pso = time.time()

    best_lh = helper.best_particle.best.likelihood

    for it in range(iterations):
        for p in particles:
            cb_particle_iteration(particle_iteration(p, helper))

    # # Uncomment the following for parallel computation
    # for it in range(iterations):
    #     print("------- Iteration %d -------" % it)

    #     start_it = time.time()

    #     pool = mp.Pool(mp.cpu_count())
    #     processes = []
    #     for p in particles:
    #         processes.append(pool.apply_async(particle_iteration, args=(p, helper), callback=cb_particle_iteration))

    #     # before starting a new iteration we wait for every process to end
    #     # for p in processes:
    #     #     p.start()
    #     #     print("Got it")
    #     pool.close()
    #     pool.join()

    #     end_it = time.time()

    #     print("Iteration %d completed in %f seconds" % (it, end_it - start_it))
    #     print("Likelihood increased by %f%%" % ((best_lh - helper.best_particle.best.likelihood) / helper.best_particle.best.likelihood * 100))

    #     if helper.best_particle.best.likelihood > best_lh:
    #         best_lh = helper.best_particle.best.likelihood

    # end_pso = time.time()

    # print("PSO completed in %f seconds" % (end_pso - start_pso))