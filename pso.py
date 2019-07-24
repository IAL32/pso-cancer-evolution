import multiprocessing as mp
import random
import sys
import time
import math

from Helper import Helper
from Node import Node
from Operation import Operation as Op
from Particle import Particle
from Tree import Tree
from datetime import datetime

from Data import Data

# global scope for multiprocessing
particles = []
helper = None
data = None

def init(nparticles, iterations, matrix, mutations, mutation_names, cells, alpha, beta, k, c1, c2, seed):
    global helper
    global particles
    global data
    helper = Helper(matrix, mutations, mutation_names, cells, alpha, beta, k, c1, c2)
    data = Data(nparticles, iterations, seed)

    pso(nparticles, iterations, matrix)
    data.helper = helper
    return data, helper

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
    # doing this for when we parallelize everything
    particles[p.number] = p

    # updating log likelihood and bests
    lh = Tree.greedy_loglikelihood(helper, tree_copy)
    tree_copy.likelihood = lh
    p.current_tree = tree_copy
    # print("Operation %d" % (tree_copy.operation.type))

    best_particle_lh = p.best.likelihood
    best_swarm_lh = helper.best_particle.best.likelihood

    if lh > best_particle_lh:
        # updating particle best
        decreased = (best_particle_lh - lh) / lh * 100
        print(mytime() + "!! %d new particle best, before: %f, now: %f, increased by %f%%" % (p.number, best_particle_lh, lh, decreased))
        data.iteration_new_particle_best[i][p.number] = lh
        p.best = tree_copy
    if lh > best_swarm_lh:
        # updating swarm best
        decreased = (best_swarm_lh - lh) / lh * 100
        print(mytime() + "!!!!! %d new swarm best, before: %f, now: %f, increased by %f%%" % (p.number, best_swarm_lh, lh, decreased))
        data.iteration_new_best[i][p.number] = lh
        helper.best_particle = p

    data.particle_iteration_times[p.number].append(data._passed_seconds(start_time, time.time()))

def particle_iteration(it, p, helper):

    start_time = time.time()
    ops = list(range(0, Op.NUMBER))
    result = -1

    op = ops.pop(random.randint(0, len(ops) - 1))

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

    tree_copy = p.current_tree.copy()
    # Higher probability of steep when no new best
    if it > 0 and data.iteration_new_particle_best[it - 1][p.number] == 0:
        p.climb_probability += 0.1
    if random.random() < math.log(p.climb_probability):
        print(mytime() + "/// %d Done operation" % p.number)
        p.climb_probability = 1.0
        result = Op.tree_operation(helper, tree_copy, op)
    else:
        print(mytime() + "/// %d Trying clade attachment" % p.number)
        best_swarm_copy = helper.best_particle.best.copy()
        best_particle_copy = p.best.copy()

        current_tree_mutations, current_tree_mn = p.current_tree.phylogeny.mutation_number(helper)
        max_clades = 2
        distance_particle, _, _, mutations_particle, mut_number_particle = p.current_tree.phylogeny.distance(helper, best_particle_copy.phylogeny)
        distance_swarm,    _, _, mutations_swarm, mut_number_swarm       = p.current_tree.phylogeny.distance(helper, best_swarm_copy.phylogeny)

        particle_clade = best_particle_copy.phylogeny.get_clade_distance(helper, max_clades, current_tree_mn, distance_particle)
        swarm_clade = best_swarm_copy.phylogeny.get_clade_distance(helper, max_clades, current_tree_mn, distance_swarm)

        if particle_clade is not None or swarm_clade is not None:
            clade_attach = None
            if distance_particle == 0 or particle_clade is None: # it is the same tree
                clade_attach = swarm_clade
            elif distance_swarm == 0 or swarm_clade is None: # it is the same tree
                clade_attach = particle_clade
            else:
                ran = random.random()
                if ran < .5:
                    clade_attach = particle_clade
                else:
                    clade_attach = swarm_clade

            clade_to_attach = tree_copy.phylogeny.get_clade_distance(helper, max_clades, current_tree_mn, max(distance_particle, distance_swarm), root=True)
            if clade_to_attach is not None:
                # print ("mut id:", clade_to_attach.mutation_id, "mut height", clade_to_attach.get_height(), "curr mt:", current_tree_mn, "dist:", max(distance_particle, distance_swarm))
                clade_attach = clade_attach.copy().detach()
                clade_to_attach.attach_clade_and_fix(helper, tree_copy, clade_attach)

            tree_copy.phylogeny.fix_for_losses(helper, tree_copy)

    return it, result, op, p, tree_copy, start_time

def particle_iteration_hill(it, p, helper):
    start_time = time.time()
    ops = list(range(0, Op.NUMBER))
    result = -1

    ran = random.random()
    if ran < .33:
        tree_copy = helper.best_particle.best.copy()
    elif ran < .66:
        tree_copy = p.best.copy()
    else:
        tree_copy = p.current_tree.copy()

    op = ops.pop(random.randint(0, len(ops) - 1))
    result = Op.tree_operation(helper, tree_copy, op)
    return it, result, op, p, tree_copy, start_time

def particle_iteration_hill_2(it, p, helper):
    start_time = time.time()
    ops = list(range(0, Op.NUMBER))
    result = -1
    # while len(ops) > 0 and result != 0:
    op = ops.pop(random.randint(0, len(ops) - 1))
    r = random.random()
    d_s, (_, highest_s) = p.current_tree.phylogeny.distance(helper, helper.best_particle.best.copy().phylogeny)
    d_p, (_, highest_p) = p.current_tree.phylogeny.distance(helper, p.best.copy().phylogeny)
    if r < .25:
        tree_copy = helper.best_particle.best.copy()
    elif r < .50:
        tree_copy = p.best.copy()
    elif r < .75:
        clade_to_attach = highest_s if d_s > d_p else highest_p
        tree_copy = Tree.random(helper.cells, helper.mutations, helper.mutation_names)
        tree_copy.phylogeny.attach_clade_and_fix(helper, tree_copy, clade_to_attach)

    else:
        tree_copy = p.current_tree.copy()
    
    result = Op.tree_operation(helper, tree_copy, op)
    return it, result, op, p, tree_copy, start_time

def particle_iteration_clades(it, p, helper):
    start_time = time.time()
    ops = list(range(0, Op.NUMBER))
    result = -1
    # while len(ops) > 0 and result != 0:
    op = ops.pop(random.randint(0, len(ops) - 1))
    r = random.random()

    best_swarm_copy = helper.best_particle.best.copy()
    best_particle_copy = p.best.copy()

    distance_particle = p.current_tree.phylogeny.distance(helper, best_particle_copy.phylogeny)
    distance_swarm = p.current_tree.phylogeny.distance(helper, best_swarm_copy.phylogeny)
    particle_clades = best_particle_copy.phylogeny.get_clades_max_nodes(max=distance_particle)
    swarm_clades = best_swarm_copy.phylogeny.get_clades_max_nodes(max=distance_swarm)

    max_clades = 2

    if distance_particle < max_clades and distance_swarm < max_clades or len(particle_clades) == 0 and len(swarm_clades) == 0:
        tree_copy = p.current_tree.copy()
        result = Op.tree_operation(helper, tree_copy, op)
    else:
        clades_attach = []

        if distance_particle == 0 or len(particle_clades) == 0:
            for i in range(max_clades):
                if len(swarm_clades) > 0:
                    choice = random.choice(swarm_clades)
                    clades_attach.append(choice)
                    swarm_clades.remove(choice)
        elif distance_particle == 0 or len(swarm_clades) == 0:
            for i in range(max_clades):
                if len(swarm_clades) > 0:
                    choice = random.choice(particle_clades)
                    clades_attach.append(choice)
                    particle_clades.remove(choice)
        else:
            for i in range(max_clades):
                ran = random.random()
                if ran < .5 and len(particle_clades) > 0:
                    # particle clade
                    choice = random.choice(particle_clades)
                elif len(swarm_clades) > 0:
                    choice = random.choice(swarm_clades)
                    clades_attach.append(choice)
                    swarm_clades.remove(choice)
        tree_copy = p.current_tree.copy()
        for cl in clades_attach:
            cl_ = cl.detach().copy()
            cta = random.choice(tree_copy.phylogeny.get_clades())
            cta.attach_clade_and_fix(helper, tree_copy, cl_)
        tree_copy.phylogeny.fix_for_losses(helper, tree_copy)

    return i, result, op, p, tree_copy, start_time

def mytime():
    return( datetime.now().strftime("[%Y/%m/%d, %H:%M:%S] - ") )

def pso(nparticles, iterations, matrix):
    global particles
    global helper
    global data
    # Particle initialization
    print(mytime() + "Particle initialization...")
    # Random position, each tree is a binary tree at the beginning
    particles = [Particle(helper.cells, helper.mutations, helper.mutation_names, n) for n in range(nparticles)]

    helper.best_particle = particles[0]
    pool = mp.Pool(mp.cpu_count())

    data.initialization_start = time.time()

    # parallelizing tree initialization
    processes = []
    for i, p in enumerate(particles):
        processes.append(pool.apply_async(init_particle, args=(i, p, helper), callback=cb_init_particle))

    pool.close()
    pool.join()

    # non-parallel initialization
    # for i, p in enumerate(particles):
    #     cb_init_particle(init_particle(i, p, helper))

    data.starting_likelihood = helper.best_particle.best.likelihood
    data.initialization_end = time.time()

    data.pso_start = time.time()

    # Uncomment the following for single core computation

    # for it in range(iterations):
    #     start_it = time.time()

    #     print("------- Iteration %d -------" % it)
    #     for p in particles:
    #         # if it == 20 and p.number == 20:
    #         #     p.current_tree.debug = True
    #         # cb_particle_iteration(particle_iteration(it, p, helper))
    #         cb_particle_iteration(particle_iteration_hill(it, p, helper)) # done!
    #         # cb_particle_iteration(particle_iteration_hill_2(it, p, helper))
    #         # cb_particle_iteration(particle_iteration_clades(it, p, helper))
    #     data.best_iteration_likelihoods.append(helper.best_particle.best.likelihood)
    #     data.iteration_times.append(data._passed_seconds(start_it, time.time()))

    # Uncomment the following for parallel computation

    for it in range(iterations):
        print(mytime() + "------- Iteration %d -------" % it)

        start_it = time.time()

        pool = mp.Pool(mp.cpu_count())
        processes = []
        for p in particles:
            processes.append(pool.apply_async(particle_iteration, args=(it, p, helper), callback=cb_particle_iteration))
            # processes.append(pool.apply_async(particle_iteration_hill, args=(it, p, helper), callback=cb_particle_iteration))
            # processes.append(pool.apply_async(particle_iteration_hill_2, args=(it, p, helper), callback=cb_particle_iteration))
            # processes.append(pool.apply_async(particle_iteration_clades, args=(it, p, helper), callback=cb_particle_iteration))

        # before starting a new iteration we wait for every process to end
        # for p in processes:
        #     p.start()
        #     print("Got it")

        pool.close()
        pool.join()

        data.best_iteration_likelihoods.append(helper.best_particle.best.likelihood)
        data.iteration_times.append(data._passed_seconds(start_it, time.time()))

    data.pso_end = time.time()
