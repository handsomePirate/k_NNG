import csv
import networkx as nx
from matplotlib import pyplot as plt

xy = dict()
edges = []
readerPoints = csv.DictReader(open('points.csv'))
readerEdges = csv.DictReader(open('edges.csv'))

i = 0
for row in readerPoints:
    xy[i] = (float(row['x']), float(row['y']))
    i += 1

for row in readerEdges:
    edges.append([float(row['e1']), float(row['e2'])])

#plt.plot(x, y)
G = nx.DiGraph()
G.add_nodes_from(xy.keys())
G.add_edges_from(edges)

for n, p in xy.items():
    G.nodes[n]['pos'] = p

nx.draw(G, xy, node_size=0.2)
plt.show()