from Node import Node

class Operation(object):

    BACK_MUTATION = 0
    DELETE_MUTATION = 1
    SWITCH_NODES = 2
    PRUNE_REGRAFT = 3

    def __init__(self, type, node_uid_1 = None, node_uid_2 = None, node_uid_3 = None):
        self.type = type
        self.node_uid_1 = node_uid_1
        self.node_uid_2 = node_uid_2
        self.node_uid_3 = node_uid_3
    
    @staticmethod
    def do(operation, helper, particle):
        if (operation.type == Operation.BACK_MUTATION):
            Operation._back_mutation(operation, helper, particle)

            particle.operations.append(Operation(operation.type, operation.node_uid_1, operation.node_uid_2, operation.node_uid_3))
        elif (operation.type == Operation.DELETE_MUTATION):
            Operation._delete_mutation(operation, helper, particle)

            particle.operations.append(Operation(operation.type, operation.node_uid_1))
        elif (operation.type == Operation.SWITCH_NODES):
            Operation._switch_nodes(operation, helper, particle)

            particle.operations.append(Operation(operation.type, operation.node_uid_1, operation.node_uid_2))
        elif (operation.type == Operation.PRUNE_REGRAFT):
            Operation._prune_regraft(operation, helper, particle)

            particle.operations.append(Operation(operation.type, operation.node_uid_1, operation.node_uid_2))

    @staticmethod
    def _back_mutation(op, helper, particle):
        node = particle.tree.find_node_by_uid(op.node_uid_1)
        candidate = particle.tree.find_node_by_uid(op.node_uid_2)
        
        node_deletion = Node(candidate.name, None, candidate.mutation_id, op.node_uid_3)

        particle.losses_list.append(node_deletion)
        particle.k_losses_list[node_deletion.mutation_id] += 1

        par = node.up
        current = node.detach()
        par.add_child(node_deletion)
        node_deletion.add_child(current)

    @staticmethod
    def _delete_mutation(op, helper, particle):
        node_delete = particle.tree.find_node_by_uid(op.node_uid_1)
        node_delete.delete_b(helper, particle)

        particle.operations.append(Operation(op.type, op.node_uid_1))

    @staticmethod
    def _switch_nodes(op, helper, particle):
        u = particle.tree.find_node_by_uid(op.node_uid_1)
        v = particle.tree.find_node_by_uid(op.node_uid_2)

        u.swap(v)

        u.fix_for_losses(helper, particle)
        v.fix_for_losses(helper, particle)

    @staticmethod
    def _prune_regraft(op, helper, particle):
        u = particle.tree.find_node_by_uid(op.node_uid_1)
        v = particle.tree.find_node_by_uid(op.node_uid_2)

        u.prune_and_reattach(v)

        u.fix_for_losses(helper, particle)

    @staticmethod
    def do_list(helper, particle, operations):
        for op in operations:
            Operation.do(op, helper, particle)