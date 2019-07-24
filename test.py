import math
from Helper import Helper
from Node import Node
import matplotlib.pyplot as plt
import networkx as nx

def plotGraph(graph,ax,title):
    pos=[(ii[1],ii[0]) for ii in graph.nodes()]
    pos_dict=dict(zip(graph.nodes(),pos))
    nx.draw(graph,pos=pos_dict,ax=ax,with_labels=True)
    ax.set_title(title)
    return

helper = Helper(None, 10, [], 10, 0, 0, 3, 0, 0)

a = Node('germline', None, -1, False)
a1 = Node('1', a, 0)
a2 = Node('2', a1, 1)
a3 = Node('3', a1, 2)
a4 = Node('4', a1, 3)
a8 = Node('8', a4, 7)
a5 = Node('5', a4, 4)
a6 = Node('6', a5, 5)
a7 = Node('7', a2, 6)
a9 = Node('9', a7, 8)

a1_b = Node('1', a6, 0, True)
a2_b = Node('2', a9, 1, True)
a10 = Node('10', a2_b, 9)

a.save()

b = Node('germline', None, -1, False)
b8 = Node('8', b, 7)
b3 = Node('3', b8, 2)
b2 = Node('2', b8, 1)
b5 = Node('5', b3, 4)
b4 = Node('4', b5, 3)
b6 = Node('6', b4, 5)
b1 = Node('1', b6, 0)
b7 = Node('7', b2, 6)
b9 = Node('9', b7, 8)

b3_b = Node('3', b1, 2, True)
b10 = Node('10', b9, 9)
b2_b = Node('8', b10, 7, True)

b.save('test2.gv')

g = a.distance(helper, b)
