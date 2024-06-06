import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import tkinter as tk
from tkinter import messagebox
import math


################### all functions required to run hexlattice.py #########################
############################ designed for sequential run ################################
################ vertex numbering, corresponding edge numbering start from  1 ###########
def show_warning(message):
    root = tk.Tk()
    root.withdraw()
    result = messagebox.askyesno("Warning", message)
    root.destroy()
    return result

def save_csv(data, file_path):
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerows(data)

def load_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            data.append([float(i) for i in row])
    return np.array(data)

def load_edge_list(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=',')
        return [tuple(map(int, row)) for row in reader]

########## hexagon vertically flat : numbering 1 from left to right & bottom to top
def generate_hexgrid(s, h, v, save_path=None):
    # Ensure h and v are even for a perfect periodic lattice
    if h % 2 != 0 or v % 2 != 0:
        raise ValueError("h and v should be even for a perfect periodic lattice.")
    
    # Check if the hexgrid file already exists
    if save_path and os.path.exists(save_path):
        hexgrid = load_csv(save_path)
        return hexgrid

    # Generate the coordinates for the hexagonal grid
    x_coords = []
    y_coords = []

    for row in range(v):
        for col in range(h):
            x = s * (col * np.sqrt(3) / 2)
            if col % 2 == 0:
                y = s * row*1.5-0.5 if row % 2 == 0 else s * row *1.5
            else:
                y = s * row*1.5-0.5 if row % 2 == 1 else s * row *1.5
            x_coords.append(x)
            y_coords.append(y)

    hexgrid = np.column_stack((x_coords, y_coords))

    # Save the hexgrid to a file
    if save_path:
        save_csv(hexgrid, save_path)

    return hexgrid

def plot_hexgrid(hexgrid, edge_list1, edge_list2, show_vertex_numbers=False):
    fig, ax = plt.subplots(figsize=(8, 8))

    # Scatter plot of the vertices
    x_values, y_values = hexgrid[:, 0], hexgrid[:, 1]
    ax.scatter(x_values, y_values)

    # Plot the first set of edges
    for edge in edge_list1:
        x1, y1 = hexgrid[edge[0]-1]
        x2, y2 = hexgrid[edge[1]-1]
        ax.plot([x1, x2], [y1, y2], color='blue')

    # Plot the second set of edges
    for edge in edge_list2:
        x1, y1 = hexgrid[edge[0]-1]
        x2, y2 = hexgrid[edge[1]-1]
        ax.plot([x1, x2], [y1, y2], color='red')

    if show_vertex_numbers:
        for i, (x, y) in enumerate(zip(x_values, y_values), start=1):
            ax.text(x, y, str(i), fontsize=10, ha='center', va='center')

    ax.set_aspect('equal', 'box')
    ax.axis('off')
    plt.show()

def euclidean_distance_matrix(points, boundary_type, box_size, save_path=None):
    num_points = len(points)
    dist_matrix = np.zeros((num_points, num_points))
    
    for i in range(num_points):
        for j in range(i + 1, num_points):
            dist = euclidean_distance(points[i], points[j], boundary_type, box_size)
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    if save_path:
        save_csv(dist_matrix, save_path)
    
    return dist_matrix

def euclidean_distance(p1, p2, boundary_type, xshift,box_size):
    delta = p1 - p2
    if boundary_type == 'x' or boundary_type == 'xy':
        delta[0] = delta[0] - box_size[0] * round(delta[0] / box_size[0])
    if boundary_type == 'y' or boundary_type == 'xy':
        delta[1] = delta[1] - box_size[1] * round(delta[1] / box_size[1])
    
    return np.sqrt((delta ** 2).sum())

def create_connectivity_list(Lx, Ly,boundary_type):
    L = Lx * Ly
    def modular_index(i,j, Lx=None, Ly=None):
        if j >= 0 and j <=L:
            return j
        else:
            return None

    connectivity_list = []
    for i in range(1,L+1):
        rownum=i//Lx
        colnum=i-rownum*Lx
        neighbors = []
        if boundary_type == 'x':
            indices = [
            ((rownum+1)%L)*Lx+(colnum-1)%L,((rownum+1)%L)*Lx+(colnum)%L,((rownum+1)%L)*Lx+(colnum+1)%L,
            ((rownum)%L)*Lx+(colnum-2)%L,((rownum)%L)*Lx+(colnum-1)%L,((rownum)%L)*Lx+(colnum+1)%L,
            ((rownum-1)%L)*Lx+(colnum-1)%L,((rownum-1)%L)*Lx+(colnum)%L,((rownum-1)%L)*Lx+(colnum+1)%L,
            ]
        elif boundary_type == 'y':
            if colnum==1:
                indices = [
                ((rownum+2)%Ly)*Lx+(colnum)%L,
                ((rownum+1)%Ly)*Lx+(colnum)%L,((rownum+1)%Lx)*Lx+(colnum+1)%L,
                ((rownum)%Ly)*Lx+(colnum+1)%L,
                ((rownum-1)%Ly)*Lx+(colnum)%L,((rownum-1)%Ly)*Lx+(colnum+1)%L,
                ((rownum-2)%Ly)*Lx+(colnum)%L,((rownum-2)%Ly)*Lx+(colnum+1)%L
                ]
            elif colnum==0:
                indices = [
                ((rownum+2)%Ly)*Lx+(colnum)%L,
                ((rownum+1)%Ly)*Lx+(colnum)%L,((rownum+1)%Lx)*Lx+(colnum-1)%L,
                ((rownum)%Ly)*Lx+(colnum-1)%L,
                ((rownum-1)%Ly)*Lx+(colnum)%L,((rownum-1)%Ly)*Lx+(colnum-1)%L,
                ((rownum-2)%Ly)*Lx+(colnum)%L,((rownum-2)%Ly)*Lx+(colnum-1)%L
                ]
            else:
                indices = [
                ((rownum+2)%Ly)*Lx+(colnum)%L,
                ((rownum+1)%Ly)*Lx+(colnum-2)%L,((rownum+1)%Lx)*Lx+(colnum-1)%L,((rownum+1)%Lx)*Lx+(colnum)%L,((rownum+1)%Lx)*Lx+(colnum+1)%L,((rownum+1)%Lx)*Lx+(colnum+2)%L,
                ((rownum)%Ly)*Lx+(colnum-2)%L,((rownum)%Lx)*Lx+(colnum-1)%L,((rownum)%Lx)*Lx+(colnum+1)%L,((rownum)%Lx)*Lx+(colnum+2)%L,
                ((rownum-1)%Ly)*Lx+(colnum-2)%L,((rownum-1)%Lx)*Lx+(colnum-1)%L,((rownum-1)%Lx)*Lx+(colnum)%L,((rownum-1)%Lx)*Lx+(colnum+1)%L,((rownum-1)%Lx)*Lx+(colnum+2)%L,
                ((rownum-2)%Ly)*Lx+(colnum)%L
                ]
        elif boundary_type == 'xy':
            indices = [
            ((rownum+1)%L)*Lx+(colnum-1)%L,((rownum+1)%L)*Lx+(colnum)%L,((rownum+1)%L)*Lx+(colnum+1)%L,
            ((rownum)%L)*Lx+(colnum-2)%L,((rownum)%L)*Lx+(colnum-1)%L,((rownum)%L)*Lx+(colnum+1)%L,
            ((rownum-1)%L)*Lx+(colnum-1)%L,((rownum-1)%L)*Lx+(colnum)%L,((rownum-1)%L)*Lx+(colnum+1)%L,                
            ((rownum+2)%Ly)*Lx+(colnum)%L,((rownum+1)%Ly)*Lx+(colnum-2)%L,((rownum+1)%Lx)*Lx+(colnum-1)%L,
            ((rownum+1)%Lx)*Lx+(colnum)%L,((rownum+1)%Lx)*Lx+(colnum+1)%L,((rownum+1)%Lx)*Lx+(colnum+2)%L,
            ((rownum)%Ly)*Lx+(colnum-2)%L,((rownum)%Lx)*Lx+(colnum-1)%L,((rownum)%Lx)*Lx+(colnum+1)%L,((rownum)%Lx)*Lx+(colnum+2)%L,
            ((rownum-1)%Ly)*Lx+(colnum-2)%L,((rownum-1)%Lx)*Lx+(colnum-1)%L,((rownum-1)%Lx)*Lx+(colnum)%L,
            ((rownum-1)%Lx)*Lx+(colnum+1)%L,((rownum-1)%Lx)*Lx+(colnum+2)%L,
            ((rownum-2)%Ly)*Lx+(colnum)%L

            ]
        else:
            indices = [
            ((rownum+2)%L)*Lx+(colnum)%L,
            ((rownum+1)%L)*Lx+(colnum-1)%L,((rownum+1)%L)*Lx+(colnum)%L,((rownum+1)%L)*Lx+(colnum+1)%L,
            ((rownum)%L)*Lx+(colnum-2)%L,((rownum)%L)*Lx+(colnum-1)%L,((rownum)%L)*Lx+(colnum+1)%L,
            ((rownum-1)%L)*Lx+(colnum-1)%L,((rownum-1)%L)*Lx+(colnum)%L,((rownum-1)%L)*Lx+(colnum+1)%L,
            ((rownum-2)%L)*Lx+(colnum)%L
            ]

        for j in indices:
            mod_index = modular_index(i,j, Lx , Ly)
            if mod_index is not None and mod_index!=i:
                neighbors.append(mod_index)
        
        connectivity_list.append((i, neighbors))
    
    return connectivity_list

def generate_and_save_edge_list(points, connectivity_list, boundary_type,xshift, box_size, threshold, filename):
    """
    Generate edge list for lattice points considering periodic boundary conditions and save it to a CSV file on the fly.

    Parameters:
        points (ndarray): Array of lattice points.
        connectivity_list (list): Connectivity list of lattice points.
        boundary_type (tuple): Boundary type for each axis ('x', 'y', 'z', or 'open').
        box_size (ndarray): Size of the periodic box along each axis.
        threshold (float): Threshold for considering an edge.
        filename (str): Name of the CSV file to save.
    """
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        #writer.writerow(['vertex1_id', 'vertex2_id', 'distance'])
        for vertex, neighbors in connectivity_list:
            vertex = int(vertex)  # Convert vertex to integer
            for neighbor in neighbors:
                neighbor = int(neighbor)  # Convert neighbor to integer
                distance = euclidean_distance(points[vertex - 1], points[neighbor - 1], boundary_type,xshift, box_size)
                if distance <= threshold and distance > 0:
                    writer.writerow([vertex, neighbor, distance])

def get_nearest_neighbors(df, tolerance=1e-6, fraction1=0.5, fraction2=0.7):    
    # Find the two smallest non-zero distances
    smallest_distances = np.round(df[:,2], 6)
    unique_distances = np.unique(smallest_distances)
    
    if len(unique_distances) < 2:
        raise ValueError("Insufficient data for two smallest distances.")
    
    # Filter edges by distance with tolerance
    edges_with_distance_1 = df[np.abs(df[:,2] - unique_distances[0]) <= tolerance]
    edges_with_distance_2 = df[np.abs(df[:,2] - unique_distances[1]) <= tolerance]
    
    # Determine the cutoff indices for edge_list1 and edge_list2 based on fractions
    cutoff_index1 = int(len(edges_with_distance_1))
    cutoff_index2 = int(len(edges_with_distance_2))
    
    # Create edge_list1 and edge_list2 based on the distances, ensuring uniqueness
    unique_edges_1 = {frozenset((int(edge[0]), int(edge[1]))) for edge in edges_with_distance_1[:cutoff_index1, :2]}
    unique_edges_2 = {frozenset((int(edge[0]), int(edge[1]))) for edge in edges_with_distance_2[:cutoff_index2, :2]}
    
    edge_list1 = [(list(edge)[0], list(edge)[1]) for edge in unique_edges_1]
    edge_list2 = [(list(edge)[0], list(edge)[1]) for edge in unique_edges_2]

    # Skip some edges based on the given fraction
    if fraction1 < 1.0:
        np.random.shuffle(edge_list1)
        edge_list1 = edge_list1[:int(len(edge_list1) * fraction1)]
    
    # Skip some edges based on the given fraction
    if fraction2 < 1.0:
        np.random.shuffle(edge_list2)
        edge_list2 = edge_list2[:int(len(edge_list2) * fraction2)]

    return edge_list1, edge_list2

def load_edge_list(file_path):
    edge_list = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            edge_list.append((int(row[0]), int(row[1])))
    return edge_list

def generate_adjacency_matrix(edge_list, num_vertices):
    adj_matrix = np.zeros((num_vertices, num_vertices))
    for edge in edge_list:
        adj_matrix[edge[0]-1, edge[1]-1] = 1
        adj_matrix[edge[1]-1, edge[0]-1] = 1  # Undirected graph
    return adj_matrix

def create_output_files(hexgrid, edge_list, output_file):
    # Create adjacency matrix
    num_vertices = len(hexgrid)
    adj_matrix = generate_adjacency_matrix(edge_list, num_vertices)

    # Find coordination numbers
    coordination_numbers = {}
    for i in range(num_vertices):
        coordination_numbers[i] = int(np.sum(adj_matrix[i]))  # Counting ones in the row

    # Generate output data
    #output_data = [['id', 'x', 'y', '<z>']] # <z> is the coordination number (# of connected neighbors)
    output_data = []
    for i, (x, y) in enumerate(hexgrid, start=1):
        output_data.append([i, x, y, coordination_numbers[i - 1]])

    # Save output data to file
    save_csv(output_data, output_file)

def generate_adjacency_list(hexgrid, edge_list, fraction):
    adjacency_list = []
    for edge in edge_list:
        v1, v2 = edge
        x1, y1 = hexgrid[v1-1]
        x2, y2 = hexgrid[v2-1]
        nx = x2 - x1 # normal vector in x of edge for Silke's code
        ny = y2 - y1 # normal vector in y of edge for Silke's code
        rij = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        normal_x = (x2 - x1) / rij
        normal_y = (y2 - y1) / rij
        m = 1 if np.random.rand() < fraction else 0  # Assign 1 with the given fraction (frictional bonds)
        adjacency_list.append([v1, v2, nx, ny, m])
    return adjacency_list

def create_adjacency_list_file(adjacency_list, file_path):
#    header = ['id1', 'id2', 'nx', 'ny', 'm'] # m is the number of double (frictional) bonds
#    adjacency_list.insert(0, header)
    save_csv(adjacency_list, file_path)