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
    with open(file, 'r') as file1:
        coords = []

        for line in file1:
            line = line.strip('{} \n')
            (a, b) = line.split(",")
            ''' 
                x and y are expressed as latitude and longitude. These are converted with the Mercator projection (from Computer assignment 1)
                into x and y coordinates.
            '''
            coord = [(float(b)*m.pi/180), (m.log((m.tan(m.pi/4+m.pi*float(a)/360))))]
            coords.append(coord)
    return np.array(coords)


def plot_points(coords, connection, path):
    '''
    Plots a map containing points for each node, blue lines for possible connections and a red thick line representing
    the closest path between two selected nodes.
    :param coords: List of the nodes and their coordinates
    :param connection: An array with all the nodes in range of eachother
    :param path: List of shortest path between the nodes

    '''
    start_time = time.time()
    #Create lines between every connection
    line = coords[connection, :]
    #line = []
    #for j, data in enumerate(connection):
    #    line.append((coords[data[0],: ], coords[data[1], :]))


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
    print('Tid för plot_points (exl.plt.show) =', time.time() - start_time, 'sekunder')
    plt.show()


def construct_graph_connection(coord_list, radie):
    """
    sorts out which nodes are in range of each other.
    :param coord_list: the coordinates of each node
    :param radie: the radius for what is considered in rang
    :return: connection: an array with all the nodes in range of each other
             connection_distance: an array with the range between all nodes which are in range of each other
    """

    connection_distance = []
    connection = []
    for j, data in enumerate(coord_list):
        '''Calculate the relative distance of the nodes'''
        distance = np.hypot(coord_list[:,0]-data[0], coord_list[:,1]-data[1])
        '''save nodes which are in range'''
        #for i, data in enumerate(distance):
        for i in range(j+1, len(distance)):
            data = distance[i]
            if data < radie:
                connection.append([j, i])
                connection_distance.append(data)


    connection_distance = np.array(connection_distance)
    connection = np.array(connection)
    return connection, connection_distance

'''
def construct_fast_graph_connection(coord_list, radie):
    """
                sorts out which nodes are in range of eachoter.
                :param coord_list: the coordinates of each node
                :param radie: the radius for what is considered in range
                :return: sparse_graph: an sparse matrix containing all the connections and their distance
                         connections: an NumPy array containing the indices which have connections
                """
    start_time = time.time()
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
    print('Tid för construct_fast_graph_connections = ',time.time()-start_time,'sekunder')
    return sparse_graph, connections
'''
def construct_fast_graph_connection(coord_list, radie):
    """
                sorts out which nodes are in range of eachoter.
                :param coord_list: the coordinates of each node
                :param radie: the radius for what is considered in range
                :return: sparse_graph: an sparse matrix containing all the connections and their distance
                         connections: an NumPy array containing the indices which have connections
                """

    connection_distance = []
    connection = []
    coord_list_tree = scipy.spatial.cKDTree(coord_list)
    for j, data in enumerate(coord_list):
        '''save nodes which are in range'''
        connections_ckd = coord_list_tree.query_ball_point(data, radie)
        for i in connections_ckd:
            #only save upper half of the matrix
            if i > j:
                    #save the connection
                    connection.append([j, i])
                    #save the relative distance of the nodes
                    connection_distance.append(np.hypot(coord_list[i,0]-data[0], coord_list[i,1]-data[1]))

    connection_distance = np.array(connection_distance)
    connection = np.array(connection)


    return connection, connection_distance


def construct_graph(indices, distances, n):
    """
    creates a csr matrix
    :param indices: indices of the values
    :param distances: the values to placed in the matrix
    :param size of matrix (N x N)
    :return: returns an sparse array of the values
    """
    CSR_graph = scipy.sparse.csr_matrix((distances, [indices[:, 0], indices[:, 1]]), shape=(n, n))
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
        path.append(j)
        j = predecessor_matrix[j]
    path.append(i)

    #reverse it so that it goes from start node to end node instead
    path.reverse()
    return path


choice = input('Run fast version? (y/n)')
if choice == 'n':
    choice2 = input('Include plotting? (y/n)')
elif choice != 'y':
    print('Invalid choice, run again')
    sys.exit(1)
else:
    choice2 = input('Include plotting? (y/n)')


#file = 'SampleCoordinates.txt'
file = 'GermanyCities.txt'
#file ='HungaryCities.txt'

if file == 'SampleCoordinates.txt':
    start_node = 0
    end_node = 5
    radie = 0.08
elif file == 'GermanyCities.txt':
    start_node = 1573
    end_node = 10584
    radie = 0.0025
else:
    start_node = 311
    end_node = 702
    radie = 0.005



tot_time = time.time()
coords = read_coordinate_file(file)
print('read_coordinate_file tid = ',time.time()-tot_time,'sekunder')
if choice == 'y':
    start_time = time.time()
    connection, connection_distance = construct_fast_graph_connection(coords, radie)
    print('construct_fast_graph_connection tid = ',time.time()- start_time,'sekunder')

else:
    start_time = time.time()
    connection, connection_distance = construct_graph_connection(coords, radie)
    print('construct_graph_connection tid = ',time.time()- start_time,'sekunder')

N = len(coords)
start_time = time.time()
csr = construct_graph(connection, connection_distance, N)
print('construct_graph tid = ',time.time()- start_time,'sekunder')
sexosjutid = time.time()
# Djikstra computes the shortest path between nodes, with no indices given it will compute the path between every
#node but if given indices it will only compute the paths to the nodes given by the indices. To save computational
#effort, the paths are only calculated to the start node.
#return_predecessors gives, if true, an NxN array with the preceddors of the node, meaning if you go to the column
#specified at the column you are currently at and do so until you've reached the the column where column=row you've
#taken the shortest path from the node represented by the starting column to the node represented by the row
#Directed =  false meanse that we tell the program that if it is possible to go from node a-b then it is also possible
#to go from b-a and it will be the same distance.
min_distances, predexessor = dijkstra(csr, return_predecessors=True, indices=start_node, directed=False)
path = compute_path(predexessor, start_node, end_node)
print('Tid för 6+7 = ', time.time()-sexosjutid,'sekunder')
if choice2 == 'y':
    plot_points(coords, connection,path)
print('Totaltid', time.time() - tot_time, 'sekunder')

print('Distance = ', min_distances[end_node])
print('best path = ', path)

