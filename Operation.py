from Node import Node

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