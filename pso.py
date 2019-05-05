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
    if (particle.last_tree().likelihood > helper.best_particle.best.likelihood):
        helper.best_particle = particle

def init_particle(i, p, helper):
    lh = Tree.greedy_loglikelihood(helper, p.last_tree())
    p.last_tree().likelihood = lh
    return i, p

def cb_particle_iteration(r):
    result, i, p, tree_copy = r
    # updating log likelihood and bests
    particles[i] = p

    if result == 0:
        lh = Tree.greedy_loglikelihood(helper, tree_copy)
        # print("Operation %d" % (tree_copy.operation.type))
        p.operations[tree_copy.operation.type] += 1
        if lh > p.best.likelihood:
            # updating particle best
            decreased = (lh - p.best.likelihood) / p.best.likelihood * 100
            print("- !! %d new particle best, before: %d, now: %d, decreased by %f%%" % (i, p.best.likelihood, lh, decreased))
            p.best = tree_copy
        if lh > helper.best_particle.best.likelihood:
            # updating swarm best
            decreased = (lh - helper.best_particle.best.likelihood) / helper.best_particle.best.likelihood * 100
            print("- !!!!! %d new swarm best, before: %d, now: %d, decreased by %f%%" % (i, helper.best_particle.best.likelihood, lh, decreased))
            helper.best_particle = p
    else:
        lh = p.last_tree().likelihood

    p.trees.append(tree_copy)
    tree_copy.likelihood = lh

def particle_iteration(i, p, helper):
    # print ("Particle n. %d" % i)
    # print ("- loglh: %d" % p.best.likelihood)
    ops = list(range(0, Op.NUMBER))
    result = -1

    op = ops.pop(random.randint(0, len(ops) - 1))
    # distance returns the distance between two trees,
    # as d(T1, T2) = n. mutations - max matching
    # and the edge connecting the two highest matching clades

    # TODO: check if the second parameter correspond to a T2 clade
    distances = []
    d_s, (_, highest_s) = p.last_tree().phylogeny.distance(helper, helper.best_particle.best.phylogeny)
    distances.append((d_s, highest_s.copy()))
    d_p, (_, highest_p) = p.last_tree().phylogeny.distance(helper, p.best.phylogeny)
    distances.append((d_p, highest_p.copy()))
    lowest_distance = (d_p, highest_p) if d_p < d_s else (d_s, highest_s)

    for p2 in particles:
        distance, (_, highest) = p.last_tree().phylogeny.distance(helper, p2.last_tree().phylogeny)
        couple = (distance, highest.copy())
        distances.append(couple)
        if distance < lowest_distance[0]:
            lowest_distance = couple

    highest_d = lowest_distance[1]
    # global > particle's > nearest
    tree_copy = Tree.random(helper.cells, helper.mutations, helper.mutation_names)
    tree_copy.phylogeny.attach_clade_and_fix(helper, tree_copy, highest_d)

    result = Op.tree_operation(helper, tree_copy, op)

    return result, i, p, tree_copy

def pso(nparticles, iterations, matrix):
    global particles
    global helper
    # Particle initialization
    print("Particle initialization...")
    # Random position, each tree is a binary tree at the beginning
    particles = [Particle(helper.cells, helper.mutations, helper.mutation_names) for n in range(nparticles)]

    helper.best_particle = particles[0]
    pool = mp.Pool(mp.cpu_count())

    start = time.time()
    processes = []
    for i, p in enumerate(particles):
        # p.last_tree().phylogeny.save("trees/tree_%d.gv" % i)
        processes.append(pool.apply_async(init_particle, args=(i, p, helper), callback=cb_init_particle))

    for p in processes:
        p.get()
    end = time.time()

    print("Initialization of %d particles: %f seconds" % (len(particles), end - start))

    pool.close()
    pool.join()

    # exit(1)

    # start = time.time()
    # for p in particles:
    #     lh = greedy_tree_loglikelihood(helper, p.last_tree())
    #     p.last_tree().likelihood = lh
    #     if (p.last_tree().likelihood > helper.best_particle.best.likelihood):
    #         helper.best = p
    # end = time.time()
    # print("Standard time for initialization: %f" % (start - end))

    print(helper.best_particle.best.likelihood)

    start_pso = time.time()
    best_lh = helper.best_particle.best.likelihood
    for it in range(iterations):
        start_it = time.time()
        print("------- Iteration %d -------" % it)
        pool = mp.Pool(mp.cpu_count())
        processes = []
        for i, p in enumerate(particles):
            processes.append(pool.apply_async(particle_iteration, args=(i, p, helper), callback=cb_particle_iteration))

        for p in processes:
            p.get()
        pool.close()
        pool.join()
        end_it = time.time()
        print("Iteration %d completed in %f seconds" % (it, end_it - start_it))
        print("Likelihood decreased by %f%%" % ((helper.best_particle.best.likelihood - best_lh) / best_lh * 100))

        if helper.best_particle.best.likelihood > best_lh:
            best_lh = helper.best_particle.best.likelihood

    end_pso = time.time()

    print("PSO completed in %f seconds" % (end_pso - start_pso))