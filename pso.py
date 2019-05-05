import copy
import math
import random as r
import sys
import time

from graphviz import Source

from Helper import Helper
from Node import Node, rid
from Operation import Operation as Op
import multiprocessing as mp
from Particle import Particle
from Tree import Tree

# r.seed(1)

particles = []
helper = None

def init(nparticles, iterations, matrix, mutations, mutation_names, cells, alpha, beta, k):
    global helper
    helper = Helper(matrix, mutations, mutation_names, cells, alpha, beta, k)

    pso(nparticles, iterations, matrix)


def cb_init_particle(result):
    i, particle = result
    particles[i] = particle
    if (particle.last_tree().likelihood > helper.best_particle.best.likelihood):
        helper.best_particle = particle

def init_particle(i, p, helper):
    lh = greedy_tree_loglikelihood(helper, p.last_tree())
    p.last_tree().likelihood = lh
    return i, p

def cb_particle_iteration(r):
    result, i, p, tree_copy = r
    # updating log likelihood and bests
    particles[i] = p
    if result == 0:
        p.operations[tree_copy.operation.type] += 1
        # print("Operation %d" % (tree_copy.operation.type))
        lh = greedy_tree_loglikelihood(helper, tree_copy)
        tree_copy.likelihood = lh
        p.trees.append(tree_copy)

        if lh > p.best.likelihood:
            # updating particle best
            decreased = (lh - p.best.likelihood) / p.best.likelihood * 100
            print("- !! Found new particle best, previous: %d, now: %d, decreased by %f%%" % (p.best.likelihood, lh, decreased))
            p.best = tree_copy
        if lh > helper.best_particle.best.likelihood:
            # updating swarm best
            decreased = (lh - helper.best_particle.best.likelihood) / helper.best_particle.best.likelihood * 100
            print("- !!!!! Found new swarm best, previous: %d, now: %d, decreased by %f%%" % (helper.best_particle.best.likelihood, lh, decreased))
            helper.best_particle = p
    

def particle_iteration(i, p, helper):
    # print ("Particle n. %d" % i)
    # print ("- loglh: %d" % p.best.likelihood)
    ops = list(range(0, Op.NUMBER))
    result = -1

    # keep trying an operation until we find a valid one
    while len(ops) > 0 and result != 0:
        # choose a random operation
        op = ops.pop(r.randint(0, len(ops) - 1))
        d_s, (_, highest_s) = p.last_tree().phylogeny.distance(helper, helper.best_particle.best.copy().phylogeny)
        d_p, (_, highest_p) = p.last_tree().phylogeny.distance(helper, p.best.copy().phylogeny)
        ran = r.random()
        if ran < .25:
            tree_copy = Tree.random(helper.cells, helper.mutations, helper.mutation_names)
            tree_copy.phylogeny.attach_clade_and_fix(helper, tree_copy, highest_s)
        elif ran < .50:
            tree_copy = Tree.random(helper.cells, helper.mutations, helper.mutation_names)
            tree_copy.phylogeny.attach_clade_and_fix(helper, tree_copy, highest_p)
        elif ran < .75:
            clade_to_attach = highest_s if d_s > d_p else highest_p
            # reducing the clade height to increase variety
            if len(clade_to_attach.children) > 0:
                clade_to_attach = clade_to_attach.get_random_node()

            tree_copy = Tree.random(helper.cells, helper.mutations, helper.mutation_names)
            tree_copy.phylogeny.attach_clade_and_fix(helper, tree_copy, clade_to_attach)
        else:
            clade_to_attach = r.choice(p.last_tree().phylogeny.get_clades())
            tree_copy = Tree.random(helper.cells, helper.mutations, helper.mutation_names)
            tree_copy.phylogeny.attach_clade_and_fix(helper, tree_copy, highest_s)
        result = tree_operation(helper, tree_copy, op)
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

    print("Multiprocessing time for initialization: %f seconds" % (end - start))

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

    for it in range(iterations):
        print("------- Iteration %d -------" % it)
        pool = mp.Pool(mp.cpu_count())
        processes = []
        for i, p in enumerate(particles):
            processes.append(pool.apply_async(particle_iteration, args=(i, p, helper), callback=cb_particle_iteration))

        for p in processes:
            p.get()
        pool.close()
        pool.join()

    for i, p in enumerate(particles):
        print ("Particle n. %d" % i)
        print ("- loglh: %d" % p.best.likelihood)
        # p.best.phylogeny.save("trees/tree_last_%d.gv" % i)
    
    print ("Overall best particle: lh: %f" % helper.best_particle.best.likelihood)
    p.best.phylogeny.save("trees/best.gv")

def tree_operation(helper, tree, operation):
    if operation == Op.BACK_MUTATION:
        # back-mutation
        return add_back_mutation(helper, tree)
    elif operation == Op.DELETE_MUTATION:
        # delete random mutation
        return mutation_delete(helper, tree)
    elif operation == Op.SWITCH_NODES:
        # switch random nodes
        return switch_nodes(helper, tree)
    elif operation == Op.PRUNE_REGRAFT:
        # prune-regraft two random nodes
        return prune_regraft(helper, tree)
    else:
        raise SystemError("Something has happened while chosing an operation")

def add_back_mutation(helper, tree):

    max_losses = helper.k
    # gets a list of all the nodes from cache
    cached_nodes = tree.phylogeny.get_cached_content()
    keys = list(cached_nodes.keys())
    # select a random node
    node = r.choice(keys)

    # node has no parent, hence cannot add a back mutation
    if (node.up == None or node.up.up == None):
        return 1
    # if losses list has reached its maximum, then we can't procede
    if (len(tree.losses_list) >= max_losses):
        return 1
    
    # select our candidates amongst the ancestors
    candidates = [p for p in node.iter_ancestors() if (p.loss == False)]

    # selecting one random ancestor
    candidate = r.choice(candidates)

    # if the ancestor is the root, we cannot procede
    if (candidate.mutation_id == -1):
        return 1
    # Ensuring we have no more than k mutations per mutation type
    if (tree.k_losses_list[candidate.mutation_id] >= helper.k):
        return 1
    # If the mutation is already lost in the current tree, no way to remove it again
    if (node.is_mutation_already_lost(candidate.mutation_id)):
        return 1
    #
    node_deletion = Node(candidate.name, None, candidate.mutation_id, rid(), True)

    tree.losses_list.append(node_deletion)
    tree.k_losses_list[node_deletion.mutation_id] += 1

    # saving parent before detaching
    par = node.up
    current = node.detach()
    par.add_child(node_deletion)
    node_deletion.add_child(current)
    current.fix_for_losses(helper, tree)
    tree.operation = Op(Op.BACK_MUTATION, node_name_1=candidate.name, node_name_2=node_deletion.name)
    return 0

def mutation_delete(helper, tree):
    if (len(tree.losses_list) == 0):
        return 1

    node_delete = r.choice(tree.losses_list)
    tree.operation = Op(Op.DELETE_MUTATION, node_name_1=node_delete.name)
    node_delete.delete_b(helper, tree)
    return 0

def switch_nodes(helper, tree):
    cached_nodes = tree.phylogeny.get_cached_content()
    keys = list(cached_nodes.keys())

    u = None
    while (u == None or u.up == None or u.loss):
        u = r.choice(keys)
        keys.remove(u)
    v = None
    keys = list(cached_nodes.keys())
    while (v == None or v.up == None or v.loss or u.name == v.name):
        v = r.choice(keys)
        keys.remove(v)

    tree.operation = Op(Op.SWITCH_NODES, node_name_1=u.name, node_name_2=v.name)

    u.swap(v)
    u.fix_for_losses(helper, tree)
    v.fix_for_losses(helper, tree)
    return 0

def prune_regraft(helper, tree):
    nodes_list = tree.phylogeny.get_cached_content()

    prune_res = -1
    while prune_res != 0:
        keys = list(nodes_list.keys())
        u = None
        while (u == None or u.up == None or u.loss):
            u = r.choice(keys)
            keys.remove(u)
        v = None
        keys = list(nodes_list.keys())
        while (v == None or v.up == None or v.loss):
            v = r.choice(keys)
            keys.remove(v)
        prune_res = u.prune_and_reattach(v)

    tree.operation = Op(Op.PRUNE_REGRAFT, node_name_1=u.name, node_name_2=v.name)
    u.fix_for_losses(helper, tree)

    return 0

def prob(I, E, genotypes, helper, particle):
    p = 0
    if I == 0:
        if E == 0:
            p = 1 - helper.beta
        elif E == 1:
            p = helper.alpha
        else:
            raise SystemError("Unknown value for E: %d" % E)
    elif I == 1:
        if E == 0:
            p = helper.beta
        elif E == 1:
            p = 1 - helper.alpha
        else:
            raise SystemError("Unknown value for E: %d" % E)
    elif I == 2:
        p = 1
    else:
        raise SystemError("Unknown value for I: %d" % I)
    return p

def accept(currentIteration, iterations):
    return r.random() < (currentIteration / iterations)

def greedy_tree_loglikelihood(helper, tree):
    "Gets maximum likelihood of a tree"
    nodes_list = tree.phylogeny.get_cached_content()
    node_genotypes = [
        [0 for j in range(helper.mutations)]
        for i in range(len(nodes_list))
    ]

    for i, n in enumerate(nodes_list):
        n.get_genotype_profile(node_genotypes[i])

    # No need to check for empty nodes, as automatically remove them

    # for i in range(node_count):
    #     if (i in range(helper.mutations)):
    #         node = nodes_list[i]
    #         node.get_genotype_profile(node_genotypes[i])
    #     else:
    #         for j in range(helper.mutations):
    #             node_genotypes[i][j] = 3

    maximum_likelihood = 0

    for i in range(helper.cells):
        best_sigma = -1
        best_lh = float("-inf")

        for n in range(len(nodes_list)):
            lh = 0
            for j in range(helper.mutations):
                p = prob(helper.matrix[i][j], node_genotypes[n][j], node_genotypes, helper, tree)
                lh += math.log(p)

            if lh > best_lh:
                best_sigma = n
                best_lh = lh
        tree.best_sigma[i] = best_sigma
        maximum_likelihood += best_lh

    return maximum_likelihood
