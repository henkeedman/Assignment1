import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import math as m
import time
import scipy.sparse
import scipy.spatial

#file='SampleCoordinates.txt'
#file = 'GermanyCities.txt'
file ='HungaryCities.txt'
if file == 'SampleCoordinates.txt':
    start_node=0
    end_node=5
    radius=0.08
elif (file == 'GermanyCities.txt'):
    start_node=10584
    end_node=1573
    radius=0.0025
else:
    start_node=311
    end_node=702
    radius=0.005

def read_coordinate_file(file):
    """             Function to read coordinatefiles
                    :param file: file to be read
                    :return: List of the nodes and their coordinates
                    """

    file1 = open(file, 'r')
    x_coords = []
    y_coords = []

    for line in file1:
        #line.readline(file1)
        line = line.replace('{', '')
        line = line.replace('}', '')
        (x, y) = line.split(",")
        y_coords.extend([m.log((m.tan(m.pi/4+m.pi*float(y)/360)))])
        x_coords.extend([float(x)*m.pi/180])

    file1.close()
    coords = np.array([y_coords, x_coords])

    return coords


def plot_points(coords, connection, path):
    line = []
    mainline=[]
    start=time.time()
    for j in range(len(connection)):
        (start, stop) = (connection[j, :])
        line.append((coords[:, int(start)], coords[:, int(stop)]))
    print(time.time()-start)
    for j in range(len(path)-1):
        start = path[j]
        stop = path[j+1]
        mainline.append((coords[:, int(start)], coords[:, int(stop)]))
    line_segments = LineCollection(line)
    mainline_segments=LineCollection(mainline, linewidths=10, colors='r')
    fig = plt.figure(figsize=(10, 15))
    plt.plot(coords[0], coords[1], 'ro', markersize=1)
    ax = fig.gca()
    ax.add_collection(line_segments)
    ax.add_collection(mainline_segments)
    ax.set_xlim((min(coords[0])), max(coords[0]))
    ax.set_ylim((min(coords[1])), max(coords[1]))
    plt.show()


def construct_graph_connection(coord_list, radius):
    """
            sorts out which nodes are in range of eachoter object.
            :param coord_list: the coordinates of each node
            :param radius: the radius for what is considered in rang
            :return: Connections; an array with all the nodes in range of eachother
                     Connection_distance; ann array with the range between all nodes which are in range of eachother
            """
    dummy = int(coord_list.size/2)
    coord_list_temp = coord_list
    connection_distance = np.array([])
    connection = np.array([])
    for j in range(dummy):
        x = coord_list[0, j]
        y = coord_list[1, j]
        '''Calculate the relative distance of the nodes'''
        coord_list_temp[0] = coord_list[0]-x
        coord_list_temp[1] = coord_list[1]-y
        distance = np.hypot(coord_list_temp[0], coord_list_temp[1])
        '''remove nodes which isn't in range'''
        for i in range(j, dummy):
            if distance[i] < radius:
                temp=[j, i], [i, j]
                connection = np.append(connection, temp)
                distance_temp=distance[i], distance[i]
                connection_distance = np.append(connection_distance, distance_temp)
    connection = connection.reshape(connection_distance.size,2)
    return connection, connection_distance



def construct_fast_graph_connection(coord_list, radius):
    """
                sorts out which nodes are in range of eachoter object.
                :param coord_list: the coordinates of each node
                :param radius: the radius for what is considered in rang
                :return: k: an array with all the nodes in range of eachother and
                 the range between all nodes which are in range of eachother
                """

    coord_list = np.transpose(coord_list)
    fortheloveofgod = scipy.spatial.cKDTree(coord_list)
    coord_lista = scipy.spatial.cKDTree(coord_list)
    k = scipy.spatial.cKDTree.sparse_distance_matrix(coord_lista, fortheloveofgod, radius, p=2.)
    return k


def construct_graph(indices, distances, N):
    """
                creates a csr matrix
                :param indices: indices of the values
                :param distances: the values to placed in the matrix
                :param size of matrix (N x N)
                :return: returns an sparse array of the values"""

    CSR_graph=scipy.sparse.csr_matrix((distances, np.transpose(indices)), shape=(N, N))
    return CSR_graph

def compute_path(predecessor_matrix, start_node, end_node):
    """
                Computes shortest path between two nodes
                :param start_node: node to go from
                :param end_node: Node to go to
                :param predecessor_matrix: predecessormatrix of the nodes
                :return: list of shortest path between the nodes
                """

    i=start_node
    j=end_node
    path = []
    while j != i :
        path=[j]+path
        j=predecessor_matrix[i,j]
    path=[i]+path
    return path

start=time.time()
coords = read_coordinate_file(file)
print(time.time()-start)
start=time.time()
(connection, connection_distance) = construct_graph_connection(coords, radius)
print(time.time()-start)
start=time.time()
#k = construct_fast_graph_connection(coords, radius)
N = coords.size/2
csr = construct_graph(connection,connection_distance, coords.size/2 )
#csr = scipy.sparse.csr_matrix(k ,shape =(N, N))
start=time.time()
print(time.time()-start)
start=time.time()
min_distances, predexessor= scipy.sparse.csgraph.dijkstra(csr,return_predecessors=True)
print(time.time()-start)

path = compute_path(predexessor, start_node, end_node)
print(time.time()-start)

plot_points(coords,connection,path)
print(time.time()-start)

print(path)