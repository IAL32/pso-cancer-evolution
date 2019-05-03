from Node import Node, rid
from graphviz import Source
import random as r
import math
import copy
# from multiprocessing import Pool as ThreadPool
from Particle import Particle
from Operation import Operation as Op
from Helper import Helper
r.seed(1)

def init(nparticles, iterations, matrix, mutations, mutation_names, cells, alpha, beta, k):
    helper = Helper(matrix, mutations, mutation_names, cells, alpha, beta, k)

    pso(nparticles, iterations, helper, matrix)

def pso(nparticles, iterations, helper, matrix):
    # Particle initialization
    print("Particle initialization...")
    # Random position, each tree is a binary tree at the beginning
    particles = [Particle(helper.cells, helper.mutations, helper.mutation_names) for n in range(nparticles)]

    helper.best_particle = particles[0]

    for i, p in enumerate(particles):
        p.last_tree.phylogeny.save("trees/tree_" + str(i) + ".gv")
        lh = greedy_tree_loglikelihood(helper, p.last_tree)
        p.last_tree.likelihood = lh
        if (lh > helper.best_particle.best.likelihood):
            helper.best = p
        print ("Particle n. %d" % i)
        print ("- loglh: %d" % lh)

    for it in range(iterations):
        print("------- Iteration %d -------" % it)
        for i, p in enumerate(particles):
            print ("Particle n. %d" % i)
            print ("- loglh: %d" % p.best.likelihood)
            op = r.random()
            tree_copy = p.last_tree.copy()
            result = tree_operation(helper, tree_copy, op)
            # successful operation
            if result == 0:
                print("- successful %d operation" % ((op * 10) % 4))
                p.trees.append(tree_copy)
                lh = greedy_tree_loglikelihood(helper, tree_copy)
                tree_copy.likelihood = lh

                # updating particle best
                if lh > p.best.likelihood:
                    p.best = p.last_tree
                    print("- !! Found new particle best")
                    # updating swarm best
                    if lh > helper.best_particle.best.likelihood:
                        helper.best_particle = p
                        print("- !!!!! Found new swarm best")

 

def tree_operation(helper, tree, operation):
    if operation < 0.25:
        # back-mutation
        back_res = add_back_mutation(helper, tree)
        if back_res == 0:
            tree.phylogeny.fix_for_losses(helper, tree)
            return 0
        return 1
    elif operation < 0.50:
        # delete random mutation
        return mutation_delete(helper, tree)
    elif operation < 0.75:
        # switch random nodes
        return switch_nodes(helper, tree)
    else:
        # prune-regraft two random nodes
        return prune_regraft(helper, tree)

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
    
    return 0

def mutation_delete(helper, tree):
    if (len(tree.losses_list) == 0):
        return 1

    node_delete = r.choice(tree.losses_list)
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
    while (v == None or v.up == None or v.loss):
        v = r.choice(keys)
        keys.remove(v)

    if u.uid == v.uid:
        return 1

    u.swap(v)

    u.fix_for_losses(helper, tree)
    v.fix_for_losses(helper, tree)

    return 0

def prune_regraft(helper, tree):
    cached_nodes = tree.phylogeny.get_cached_content()
    keys = list(cached_nodes.keys())

    prune_res = 0
    pruned_node = None
    while prune_res != 0:
        u = None
        while (u == None or u.up == None or u.loss):
            u = r.choice(keys)
            keys.remove(u)
        v = None
        keys = list(cached_nodes.keys())
        while (v == None or v.up == None or v.loss):
            v = r.choice(keys)
            keys.remove(v)
        prune_res = u.prune_and_reattach(v)
        if prune_res == 0:
            pruned_node = u
    if pruned_node is not None:
        pruned_node.fix_for_losses(helper, tree)
        return 0
    return 1

def prob(I, E, genotypes, helper, particle):
    p = 0
    if I == 0:
        if E == 0:
            p = 1 - helper.beta
        elif E == 1:
            p = helper.alpha
        else:
            raise SystemError("Unknown value for E: " + str(E))
    elif I == 1:
        if E == 0:
            p = helper.beta
        elif E == 1:
            p = 1 - helper.alpha
        else:
            raise SystemError("Unknown value for E: " + str(E))
    elif I == 2:
        p = 1
    else:
        raise SystemError("Unknown value for I")
    return p

def accept(currentIteration, iterations):
    return r.random() < (currentIteration / iterations)

def greedy_tree_loglikelihood(helper, tree):
    "Gets maximum likelihood of a tree"
    nodes_list = [n for n in tree.phylogeny.traverse()]
    node_genotypes = [
        [0 for j in range(helper.mutations)]
        for i in range(len(nodes_list))
    ]

    for i, n in enumerate(nodes_list):
        n.get_genotype_profile(node_genotypes[i])

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

        for n, node in enumerate(nodes_list):
            if node_genotypes[n][0] != 3:
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
