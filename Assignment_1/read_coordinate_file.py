import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import math as m
import time
import scipy.sparse
import scipy.spatial
from scipy.sparse.csgraph import dijkstra
import sys


def read_coordinate_file(file):
    """             Function to read coordinatefiles
                    :param file: file to be read
                    :return coords: NumPy array of the nodes and their coordinates
    """

    file1 = open(file, 'r')
    coords = []

    for line in file1:
        line = line.strip('{} \n')
        (y, x) = line.split(",")
        ''' 
            x and y are expressed as latitude and longitude. These are converted with the Mercator projection (from Computer assignment 1)
            into x and y coordinates.
        '''
        coord = [(float(x)*m.pi/180), (m.log((m.tan(m.pi/4+m.pi*float(y)/360))))]
        coords.append(coord)
    file1.close()
    print(coords)
    return np.array(coords)


def plot_points(coords, connection, path):
    '''
    Plots a map containing points for each node, blue lines for possible connections and a red thick line representing
    the closest path between two selected nodes.
    :param coords: List of the nodes and their coordinates
    :param connection: An array with all the nodes in range of eachother
    :param path: List of shortest path between the nodes

    '''
    line = []
    check = time.time()

    #Create lines between every connection
    for j, data in enumerate(connection):
        line.append((coords[data[0],: ], coords[data[1], :]))

    #Create the line showing the shortest path between two nodes
    mainline = [coords[path, :]]
    line_segments = LineCollection(line)
    mainline_segments = LineCollection(mainline, linewidths=10, colors='r')
    fig = plt.figure(figsize=(10, 15))
    plt.plot(coords[:,0], coords[:,1], 'ro', markersize=1)
    ax = fig.gca()
    ax.add_collection(line_segments)
    ax.add_collection(mainline_segments)
    ax.set_xlim((min(coords[0])), max(coords[0]))
    ax.set_ylim((min(coords[1])), max(coords[1]))
    plt.axis('Equal')
    print('plottid=', time.time() - check, 'sekunder')
    plt.show()


def construct_graph_connection(coord_list, radie):
    """
            sorts out which nodes are in range of each other object.
            :param coord_list: the coordinates of each node
            :param radie: the radius for what is considered in rang
            :return: connection: an array with all the nodes in range of eachother
                     connection_distance: ann array with the range between all nodes which are in range of eachother
            """

    connection_distance = []
    connection = []
    for j, data in enumerate(coord_list):
        x = data[0]
        y = data[1]

        '''Calculate the relative distance of the nodes'''
        distance = np.hypot(coord_list[:, 0]-x, coord_list[:, 1]-y)

        '''add nodes which are in range'''
        for i, data in enumerate(distance):
            if data < radie:
                connection.append([i, j])
                connection_distance.append(data)

    connection_distance = np.array(connection_distance)
    connection = np.array(connection)
    print(connection)
    return connection, connection_distance


def construct_fast_graph_connection(coord_list, radie):
    """
                sorts out which nodes are in range of eachoter object.
                :param coord_list: the coordinates of each node
                :param radie: the radius for what is considered in range
                :return: csr: an sparse matrix containing all the connections and their distance
                         connections: an NumPy array containing the indices which have connections
                """

    coord_list_tree = scipy.spatial.cKDTree(coord_list)
    sparse_graph = scipy.spatial.cKDTree.sparse_distance_matrix(coord_list_tree, coord_list_tree, radie, p=2.)
    connections_ckd = coord_list_tree.query_ball_tree(coord_list_tree, radie)
    connections = []

    #Changes the matrix format to match the one generated in the slower version
    # (to allow for it to be used the same way)
    for j, data in enumerate(connections_ckd):
        for item in data:
            connections.append([j, item])
    connections = np.array(connections)
    connections = connections.reshape(len(connections), 2)

    return sparse_graph, connections


def construct_graph(indices, distances, N):
    """
                creates a csr matrix
                :param indices: indices of the values
                :param distances: the values to placed in the matrix
                :param size of matrix (N x N)
                :return: returns an sparse array of the values"""

    CSR_graph=scipy.sparse.csr_matrix((distances, indices), shape=(N, N))
    return CSR_graph


def compute_path(predecessor_matrix, start_node, end_node):
    """
                Computes shortest path between two nodes
                :param start_node: node to go from
                :param end_node: Node to go to
                :param predecessor_matrix: predecessormatrix of the nodes, from Dijikstra
                :return: list of shortest path between the nodes
                """

    i = start_node
    j = end_node
    path = []

    #Go through the predecessor matrix to save the data in a list
    while j != i:
        path = [j]+path

        j = predecessor_matrix[0, j]

    path = [i]+path
    return path


choice = input('Run fast version? (y/n)')
if choice == 'n':
    choice2 = input('Include plotting? (y/n)')
elif choice != 'y':
    print('Invalid choice, run again')
    sys.exit(1)
else:
    choice2 = input('Include plotting? (y/n)')

file = 'SampleCoordinates.txt'
#file = 'GermanyCities.txt'
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


start_time = time.time()
coords = read_coordinate_file(file)


if choice == 'y':
    snabb_tid = time.time()
    csr, connection = construct_fast_graph_connection(coords, radie)
    print('Tid för CkdTree = ',  time.time()-snabb_tid,'sekunder')
    min_distances, predexessor = dijkstra(csr, return_predecessors=True, indices=[start_node])
    path = compute_path(predexessor, start_node, end_node)
    print(coords[path, 0],coords[path, 1])
    if choice2 == 'y':
        plot_points(coords,connection,path)
    print('Totaltid', time.time()-start_time,'sekunder')

else:
    connection, connection_distance = construct_graph_connection(coords, radie)
    N = coords.size/2
    csr = construct_graph(connection, connection_distance, N)
    min_distances, predexessor = scipy.sparse.csgraph.dijkstra(csr, return_predecessors=True, indices=[start_node])
    path = compute_path(predexessor, start_node, end_node)
    print('Totaltid exl. plot', time.time() - start_time, 'sekunder')
    if choice2 == 'y':
        plot_points(coords, connection,path)

print('Distance = ', min_distances[0, end_node])
print('best path = ', path)

