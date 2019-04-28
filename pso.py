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

    # initialization
    particles = [
        Particle.random(helper)
        for n in range(nparticles)
    ]

    print("Particle initialization...")
    for j, p in enumerate(particles):
        p.tree.save("trees/tree_" + str(j) + ".gv")
        p.best_likelihood = greedy_tree_loglikelihood(helper, p)
        if (p.best_likelihood > helper.best_likelihood):
            helper.best_likelihood = p.best_likelihood
            helper.best_likelihood_particle = p
        print ("Particle n. " + str(j))
        print ("- loglh: " + str(p.best_likelihood))
    
    for i in range(iterations):
        for j, p in enumerate(particles):
            op = r.random()
            tmp_particle = p.copy()
            result = particle_operation(helper, tmp_particle, op)
            # successful operation
            if result == 0:
                lh = greedy_tree_loglikelihood(helper, tmp_particle)

                # updating particle best
                if lh > p.best_likelihood:
                    p.best_likelihood = lh
                    p.best_tree = tmp_particle.tree
                    p.losses_list = tmp_particle.losses_list
                    p.k_losses_list = tmp_particle.k_losses_list

                    # updating swarm best
                    if lh > helper.best_likelihood:
                        helper.best_likelihood = lh
                        helper.best_likelihood_particle = tmp_particle

                particles[j] = Particle.random(helper)

                if accept(i, iterations): # going for the global's best
                    # I take the clades up until a maximum height
                    # depending on the current iteration ratio.
                    # We get to higher clades as the time goes on
                    height = math.ceil(tmp_particle.tree.get_height() * (i / iterations))
                    clades = tmp_particle.tree.get_clades_at_height(height)
                else: # going for the global's
                    height = math.ceil(helper.best_likelihood_particle.tree.get_height() * (i / iterations))
                    clades = helper.best_likelihood_particle.tree.get_clades_at_height(height)

                # now that we have chosen wether to go for the
                # particle's best or global best, we update the particle's
                # tree regardless
                reattach_clades(helper, particles[j], clades)

    # for i in range(iterations):

def particle_operation(helper, particle, operation):
    if operation < 0.25:
        # back-mutation
        back_res = add_back_mutation(helper, particle)
        if back_res == 0:
            particle.tree.fix_for_losses(helper, particle)
            return 0
        return 1
    elif operation < 0.50:
        # delete random mutation
        return mutation_delete(helper, particle)
    elif operation < 0.75:
        # switch random nodes
        return switch_nodes(helper, particle)
    else:
        # prune-regraft two random nodes
        return prune_regraft(helper, particle)

def add_back_mutation(helper, particle):

    max_losses = helper.k
    # gets a list of all the nodes from cache
    cached_nodes = particle.tree.get_cached_content()
    keys = list(cached_nodes.keys())
    node = r.choice(keys)

    if (node.up == None or node.up.up == None):
        return 1
    if (len(particle.losses_list) >= max_losses):
        return 1
    candidates = [p for p in node.iter_ancestors() if (p.loss == False)]

    random_deletion = r.randint(0, len(candidates) - 2) + 1
    candidate = candidates[random_deletion]
    if (candidate.mutation_id == -1):
        return 1
    # Ensuring we have no more than k mutations per mutation type
    if (particle.k_losses_list[candidate.mutation_id] >= helper.k):
        return 1
    # If the mutation is already lost in the current tree, no way to remove it again
    if (node.is_mutation_already_lost(candidate.mutation_id)):
        return 1
    #
    node_deletion = Node(candidate.name, None, candidate.mutation_id, rid(), True)

    particle.losses_list.append(node_deletion)
    particle.k_losses_list[node_deletion.mutation_id] += 1

    # saving parent before detaching
    par = node.up
    current = node.detach()
    par.add_child(node_deletion)
    node_deletion.add_child(current)
    
    return 0

def mutation_delete(helper, particle):
    if (len(particle.losses_list) == 0):
        return 1

    node_delete = r.choice(particle.losses_list)
    node_delete.delete_b(helper, particle)
    return 0

def switch_nodes(helper, particle):
    cached_nodes = particle.tree.get_cached_content()
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

    u.fix_for_losses(helper, particle)
    v.fix_for_losses(helper, particle)

    return 0

def prune_regraft(helper, particle):
    cached_nodes = particle.tree.get_cached_content()
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
        pruned_node.fix_for_losses(helper, particle)
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

def greedy_tree_loglikelihood(helper, particle):
    cached_nodes = particle.tree.get_cached_content()

    node_count = len(cached_nodes)
    nodes_list = list(cached_nodes.keys())
    node_genotypes = [
        [0 for j in range(helper.mutations)]
        for i in range(node_count)
    ]

    for i in range(node_count):
        if (i in range(helper.mutations)):
            node = nodes_list[i]
            node.get_genotype_profile(node_genotypes[i])
        else:
            for j in range(helper.mutations):
                node_genotypes[i][j] = 3

    maximum_likelihood = 0

    for i in range(helper.cells):
        best_sigma = -1
        best_lh = float("-inf")

        for n in range(node_count):
            if node_genotypes[n][0] != 3:
                lh = 0
                for j in range(helper.mutations):
                    p = prob(helper.matrix[i][j], node_genotypes[n][j], node_genotypes, helper, particle)
                    lh += math.log(p)

                if lh > best_lh:
                    best_sigma = n
                    best_lh = lh
        particle.best_sigma[i] = best_sigma
        maximum_likelihood += best_lh
    
    return maximum_likelihood

def reattach_clades(helper, particle, clades):

    for cl in clades:
        clade_nodes = cl.get_cached_content()
        # remove every node already in clades
        for cn in clade_nodes:
            to_delete = next(particle.tree.iter_search_nodes(name=cn.name), None)
            if to_delete is not None:
                to_delete.delete()

    leaves = particle.tree.get_leaves()
    candidate_leaves = leaves.copy()

    # less leaves that i can re-attach to
    if len(candidate_leaves) < len(clades):
        # select more random leaves to attach, reusing them
        candidate_leaves.append(r.choices(leaves, k = len(clades) - len(candidate_leaves)))
    
    for leaf in candidate_leaves:
        for cl in clades:
            leaf.add_child(cl)

    losses_list, k_losses_list = particle.calculate_losses_list(helper)
    particle.losses_list = losses_list
    particle.k_losses_list = k_losses_list
    particle.tree.fix_for_losses(helper, particle)