from Node import Node
import random as r

def accept(currentIteration, iterations):
    return r.random() < (currentIteration / iterations)

class Operation(object):

    BACK_MUTATION = 0
    DELETE_MUTATION = 1
    SWITCH_NODES = 2
    PRUNE_REGRAFT = 3

    NUMBER = 4

    def __init__(self, type, node_name_1 = None, node_name_2 = None, node_name_3 = None):
        self.type = type
        self.node_name_1 = node_name_1
        self.node_name_2 = node_name_2
        self.node_name_3 = node_name_3

    @classmethod
    def tree_operation(cls, helper, tree, operation):
        if operation == cls.BACK_MUTATION:
            # back-mutation
            return cls.add_back_mutation(helper, tree)
        elif operation == cls.DELETE_MUTATION:
            # delete random mutation
            return cls.mutation_delete(helper, tree)
        elif operation == cls.SWITCH_NODES:
            # switch random nodes
            return cls.switch_nodes(helper, tree)
        elif operation == cls.PRUNE_REGRAFT:
            # prune-regraft two random nodes
            return cls.prune_regraft(helper, tree)
        else:
            raise SystemError("Something has happened while chosing an operation")

    @classmethod
    def add_back_mutation(cls, helper, tree):

        max_losses = helper.k
        # gets a list of all the nodes from cache
        cached_nodes = tree.phylogeny.get_cached_content()
        keys = list(cached_nodes.keys())
        # select a random node
        # node has no parent, hence cannot add a back mutation
        # keep trying until we find a suitable node
        node = r.choice(keys)
        if node.up == None or node.up.up == None:
            return 1

        # if losses list has reached its maximum, then we can't procede
        if (len(tree.losses_list) >= max_losses):
            return 1

        # select our candidates amongst the ancestors
        candidates = [p for p in node.iter_ancestors() if (p.loss == False) and (p.mutation_id != -1)]
        if len(candidates) == 0:
            return 1

        # selecting one random ancestor
        candidate = r.choice(candidates)
        # Ensuring we have no more than k mutations per mutation type
        if (tree.k_losses_list[candidate.mutation_id] >= helper.k):
            return 1
        # If the mutation is already lost in the current tree, no way to remove it again
        if (node.is_mutation_already_lost(candidate.mutation_id)):
            return 1
        #
        node_deletion = Node(candidate.name, None, candidate.mutation_id, True)

        tree.losses_list.append(node_deletion)
        tree.k_losses_list[node_deletion.mutation_id] += 1

        # saving parent before detaching
        par = node.up
        current = node.detach()
        par.add_child(node_deletion)
        node_deletion.add_child(current)
        current.fix_for_losses(helper, tree)
        tree.operation = cls(cls.BACK_MUTATION, node_name_1=candidate.name, node_name_2=node_deletion.name)
        return 0

    @classmethod
    def mutation_delete(cls, helper, tree):
        if (len(tree.losses_list) == 0):
            return 1

        node_delete = r.choice(tree.losses_list)
        tree.operation = cls(cls.DELETE_MUTATION, node_name_1=node_delete.name)
        node_delete.delete_b(helper, tree)
        return 0

    @classmethod
    def switch_nodes(cls, helper, tree):
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

        tree.operation = cls(cls.SWITCH_NODES, node_name_1=u.name, node_name_2=v.name)

        u.swap(v)
        u.fix_for_losses(helper, tree)
        v.fix_for_losses(helper, tree)
        return 0

    @classmethod
    def prune_regraft(cls, helper, tree):
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

        tree.operation = cls(cls.PRUNE_REGRAFT, node_name_1=u.name, node_name_2=v.name)
        u.fix_for_losses(helper, tree)

        return 0

    @classmethod
    def prob(cls, I, E, genotypes, helper, particle, data=None):
        p = 0
        if I == 0:
            if E == 0:
                # TODO: sigma
                if data is not None:
                    data.true_negative += 1
                p = 1 - helper.beta
            elif E == 1:
                if data is not None:
                    data.false_negatives += 1
                p = helper.alpha
            else:
                raise SystemError("Unknown value for E: %d" % E)
        elif I == 1:
            if E == 0:
                if data is not None:
                    data.false_positives += 1
                p = helper.beta
            elif E == 1:
                if data is not None:
                    data.true_positive += 1
                p = 1 - helper.alpha
            else:
                raise SystemError("Unknown value for E: %d" % E)
        elif I == 2:
            if data:
                data.missing_values += 1
            p = 1
        else:
            raise SystemError("Unknown value for I: %d" % I)
        return p
