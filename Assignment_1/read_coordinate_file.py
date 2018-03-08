import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import math as m
import time
import scipy.sparse
import scipy.spatial

#file = 'SampleCoordinates.txt'
file = 'GermanyCities.txt'
#file ='HungaryCities.txt'

if file == 'SampleCoordinates.txt':
    start_node = 0
    end_node = 5
    radie = 0.08
elif file == 'GermanyCities.txt':
    start_node = 10584
    end_node = 1573
    radie = 0.0025
else:
    start_node = 311
    end_node = 702
    radie = 0.005


def read_coordinate_file(file):
    """             Function to read coordinatefiles
                    :param file: file to be read
                    :return: List of the nodes and their coordinates
                    """

    file1 = open(file, 'r')
    x_coords = []
    y_coords = []

    for line in file1:
        line = line.strip()
        line = line.strip('{')
        line = line.strip('}')
        (y, x) = line.split(",")
        line.strip()
        ''' 
            x and y are expressed as latitude and longitude. These are converted with the Mercator projection (from Computer assignment 1)
            into x and y coordinates.
        '''

        y_coords.extend([m.log((m.tan(m.pi/4+m.pi*float(y)/360)))])
        x_coords.extend([float(x)*m.pi/180])

    file1.close()
    coords = np.array([x_coords, y_coords])
    return coords


def plot_points(coords, connection, path):
    line = []
    mainline = []
    starttime=time.time()

    for j, data in enumerate(connection):
        start, stop = data
        line.append((coords[:, start], coords[:, stop]))

    start = path[0]
    for j, data in enumerate(path[1:-1]):
        stop = data
        mainline.append((coords[:, start], coords[:, stop]))
        start = stop

    line_segments = LineCollection(line)
    mainline_segments = LineCollection(mainline, linewidths=10, colors='r')
    fig = plt.figure(figsize=(10, 15))
    plt.plot(coords[0], coords[1], 'ro', markersize=1)
    ax = fig.gca()
    ax.add_collection(line_segments)
    ax.add_collection(mainline_segments)
    ax.set_xlim((min(coords[0])), max(coords[0]))
    ax.set_ylim((min(coords[1])), max(coords[1]))
    print('plotting time = ', time.time() - starttime)
    plt.show()


def construct_graph_connection(coord_list, radius):
    """
            sorts out which nodes are in range of each other object.
            :param coord_list: the coordinates of each node
            :param radius: the radius for what is considered in rang
            :return: Connections; an array with all the nodes in range of eachother
                     Connection_distance; ann array with the range between all nodes which are in range of eachother
            """
    dummy = int(coord_list.size/2)
    coord_list_temp = coord_list
    connection_distance = []
    connection = []
    for j, data in enumerate(coord_list.transpose()):
        x = data[0]
        y = data[1]

        '''Calculate the relative distance of the nodes'''
        coord_list_temp[0] = coord_list[0]-x
        coord_list_temp[1] = coord_list[1]-y
        distance = np.hypot(coord_list_temp[0], coord_list_temp[1])

        '''remove nodes which isn't in range'''
        for i, data in enumerate(distance):
            if data < radie:
                connection.append([i, j])
                connection_distance.append(data)

    connection_distance = np.array(connection_distance)
    connection = np.array(connection)
    connection = connection.reshape(len(connection_distance), 2)
    return connection, connection_distance


def construct_fast_graph_connection(coord_list, radie):
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
    k = scipy.spatial.cKDTree.sparse_distance_matrix(coord_lista, fortheloveofgod, radie, p=2.)
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

    i = start_node
    j = end_node
    path = []
    while j != i :
        path = [j]+path
        j = predecessor_matrix[i,j]
    path = [i]+path
    return path


start = time.time()
coords = read_coordinate_file(file)
#print(time.time()-start)

grafen=time.time()
(connection, connection_distance) = construct_graph_connection(coords, radie)
#print(time.time()-start)
grafen=time.time() - grafen
print('construct_graph_time = ', grafen)

#k = construct_fast_graph_connection(coords, radius)
N = coords.size/2
csr = construct_graph(connection,connection_distance, coords.size/2 )
#csr = scipy.sparse.csr_matrix(k ,shape =(N, N))

min_distances, predexessor= scipy.sparse.csgraph.dijkstra(csr,return_predecessors=True)
#print(time.time()-start)



path = compute_path(predexessor, start_node, end_node)
#print(time.time()-start)

plot_points(coords, connection, path)
#print(time.time()-start)
print(min_distances[start_node, end_node])
print(path)