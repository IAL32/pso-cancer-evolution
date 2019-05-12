import multiprocessing as mp
import random
import sys
import time

from Helper import Helper
from Node import Node, rid
from Operation import Operation as Op
from Particle import Particle
from Tree import Tree

from Data import Data

# random.seed(1)

# global scope for multiprocessing
particles = []
helper = None
data = None

def init(nparticles, iterations, matrix, mutations, mutation_names, cells, alpha, beta, k):
    global helper
    global particles
    global data
    helper = Helper(matrix, mutations, mutation_names, cells, alpha, beta, k)
    data = Data(nparticles, iterations)

    pso(nparticles, iterations, matrix)
    data.summary(helper)

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
    i, result, op, p, tree_copy, start_time = r
    # updating log likelihood and bests
    particles[p.number] = p

    if result == 0:
        lh = Tree.greedy_loglikelihood(helper, tree_copy)
        tree_copy.likelihood = lh
        p.current_tree = tree_copy
        # print("Operation %d" % (tree_copy.operation.type))
        if lh != p.best.likelihood and lh > p.best.likelihood:
            # updating particle best
            decreased = (p.best.likelihood - lh) / lh * 100
            print("- !! %d new particle best, before: %f, now: %f, increased by %f%%" % (p.number, p.best.likelihood, lh, decreased))
            data.iteration_new_particle_best[i][p.number] = lh
            p.best = tree_copy
        if lh > helper.best_particle.best.likelihood:
            # updating swarm best
            decreased = (helper.best_particle.best.likelihood - lh) / lh * 100
            print("- !!!!! %d new swarm best, before: %f, now: %f, increased by %f%%" % (p.number, helper.best_particle.best.likelihood, lh, decreased))
            data.iteration_new_best[i][p.number] = lh
            helper.best_particle = p

    data.particle_iteration_times[p.number].append(data._passed_seconds(start_time, time.time()))

def particle_iteration(it, p, helper):

    start_time = time.time()
    ops = list(range(0, Op.NUMBER))
    result = -1

    # while len(ops) > 0 and result != 0:

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

    distance_particle = p.current_tree.phylogeny.distance(helper, best_particle_copy.phylogeny)
    distance_swarm = p.current_tree.phylogeny.distance(helper, best_swarm_copy.phylogeny)

    particle_clades = best_particle_copy.phylogeny.get_clades_max_nodes(max=distance_particle)
    swarm_clades = best_swarm_copy.phylogeny.get_clades_max_nodes(max=distance_swarm)

    max_clades = 3

    if len(particle_clades) == 0 and len(swarm_clades) == 0:
        tree_copy = p.current_tree.copy()
    else:
        clades_attach = []

        if distance_particle == 0 or len(particle_clades) == 0: # it is the same tree
            for i in range(max_clades):
                if len(swarm_clades) > 0:
                    choice = random.choice(swarm_clades)
                    clades_attach.append(choice)
                    swarm_clades.remove(choice)
        elif distance_swarm == 0 or len(swarm_clades) == 0: # it is the same tree
            for i in range(max_clades):
                if len(particle_clades) > 0:
                    choice = random.choice(particle_clades)
                    clades_attach.append(choice)
                    particle_clades.remove(choice)
        else:
            for i in range(max_clades):
                ran = random.random()
                if ran < .5 and len(particle_clades) > 0:
                    # particle clade
                    choice = random.choice(particle_clades)
                    clades_attach.append(choice)
                    particle_clades.remove(choice)
                elif ran >= .5 and len(swarm_clades) > 0:
                    # swarm clade
                    choice = random.choice(swarm_clades)
                    clades_attach.append(choice)
                    swarm_clades.remove(choice)

        tree_copy = p.current_tree.copy()

        # attaching the clades we selected before
        # and fixing every time we add a new clade
        for cl in clades_attach:
            cl_ = cl.detach().copy()
            tree_copy_clades = tree_copy.phylogeny.get_clades_at_average_level(percent=it / data.iterations)
            clade_to_attach = random.choice(tree_copy_clades)
            clade_to_attach.attach_clade_and_fix(helper, tree_copy, cl_)

        tree_copy.phylogeny.fix_for_losses(helper, tree_copy)

    result = Op.tree_operation(helper, tree_copy, op)

    return it, result, op, p, tree_copy, start_time

def pso(nparticles, iterations, matrix):
    global particles
    global helper
    global data
    # Particle initialization
    print("Particle initialization...")
    # Random position, each tree is a binary tree at the beginning
    particles = [Particle(helper.cells, helper.mutations, helper.mutation_names, n) for n in range(nparticles)]

    helper.best_particle = particles[0]
    pool = mp.Pool(mp.cpu_count())

    data.initialization_start = time.time()
    processes = []
    # parallelizing tree initialization
    for i, p in enumerate(particles):
        processes.append(pool.apply_async(init_particle, args=(i, p, helper), callback=cb_init_particle))

    for p in processes:
        p.get()

    data.initialization_end = time.time()
    data.starting_likelihood = helper.best_particle.best.likelihood

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

    data.pso_start = time.time()

    # Uncomment the following for single core computation

    for it in range(iterations):
        start_it = time.time()

        print("------- Iteration %d -------" % it)
        for p in particles:
            cb_particle_iteration(particle_iteration(it, p, helper))

        data.iteration_times.append(data._passed_seconds(start_it, time.time()))

    # Uncomment the following for parallel computation

    # for it in range(iterations):
    #     print("------- Iteration %d -------" % it)

    #     start_it = time.time()

    #     pool = mp.Pool(mp.cpu_count())
    #     processes = []
    #     for p in particles:
    #         processes.append(pool.apply_async(particle_iteration, args=(it, p, helper), callback=cb_particle_iteration))

    #     # before starting a new iteration we wait for every process to end
    #     # for p in processes:
    #     #     p.start()
    #     #     print("Got it")

    #     pool.close()
    #     pool.join()

    #     data.iteration_times.append(data._passed_seconds(start_it, time.time()))

    data.pso_end = time.time()
