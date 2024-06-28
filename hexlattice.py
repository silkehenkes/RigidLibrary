import numpy as np
import os
import pandas as pd
import csv
from hexlattice_func import *

# Example usage
s = 1  # size of the hexagon


# Lattice size
nhex = 100


# Lattice aspect ratio
aratio=.4

# Length aspect ratio
lratio=np.sqrt(3.)/2.*aratio


print("Hexagonal lattice with lattice aspect ratio ",aratio)
print(" Lx = ",nhex*np.sqrt(3)/2.," Ly = ",nhex*aratio*3./2," Lx/Ly = ",1/lratio,"\n")


# Fraction settings
fraction_edge1 = 0.85
fraction_edge2 = 0.25
fraction_dbonds = 0.5  # Adjust the fraction here

# Define boundary type
boundary_type = 'x'  # options: 'open', 'x', 'y', 'xy'

# Options to run specific parts
run_hexgrid = True
run_connectivity_list = True
plot_grid_array= True
run_edge_lists = True
run_adjacency_list = True

topdir='./KyungLattice/'

def main():
    try:
        h = nhex
        v = int(nhex*aratio)

        ######### Please check if these files are in your running directory if you want to generate new one
        hexgrid_path = topdir+f'hexgrid_{h}by{v}_{boundary_type}.csv'
        connectivity_list_path = topdir+f'connectivity_list_{h}by{v}_{boundary_type}.csv'
        edge_list_path = topdir+f'edge_list_{h}by{v}_{boundary_type}.csv'
        edge_list1_path = topdir+f'edge_list1_{h}by{v}_{boundary_type}.csv'
        edge_list2_path = topdir+f'edge_list2_{h}by{v}_{boundary_type}.csv'
        output_file = topdir+f'particle_positions_{h}by{v}_{boundary_type}.txt'
        adjacency_list_path = topdir+f'Adjacency_list_{h}by{v}_{boundary_type}.txt'

        # Generate hexgrid
        if run_hexgrid:
            if os.path.exists(hexgrid_path):
                if show_warning("Are you using previous hexgrid.csv as an input?"):
                    hexgrid = load_csv(hexgrid_path)
                else:
                    print("Terminated by user.")
                    return
            else:
                hexgrid = generate_hexgrid(s, h, v, save_path=hexgrid_path)
        else:
            hexgrid = load_csv(hexgrid_path)

        # Define box size
        box_size = [s * h* np.sqrt(3)/2, s * v *0.5+s * (v-2)] #### for hexgrid
        xshift=s * (h)* np.sqrt(3)/2*aratio
        threshold = 2  # Set the threshold for considering an edge

        # Generate the connectivity & edge list and save it
        if run_connectivity_list :
            if os.path.exists(connectivity_list_path):
                print("Overwriting connectivity_list.csv.")
                list=create_connectivity_list(h,v,boundary_type)
                save_csv(list, connectivity_list_path)
                # Save edge list to CSV file
                generate_and_save_edge_list(hexgrid,list, boundary_type, xshift,box_size, threshold, edge_list_path)
            else:
                list=create_connectivity_list(h,v,boundary_type)
                save_csv(list, connectivity_list_path)
                # Save edge list to CSV file
                generate_and_save_edge_list(hexgrid,list, boundary_type, xshift,box_size, threshold, edge_list_path)
        else:
            return
        
        # Create the edge lists and save them
        if run_edge_lists:
            if os.path.exists(edge_list1_path):
                if show_warning("Are you using previous edge_list1.csv as an input?"):
                    edge_list = load_csv(edge_list_path)
                    edge_list1, edge_list2 = get_nearest_neighbors(edge_list, tolerance=1e-6,fraction1=fraction_edge1,fraction2=fraction_edge2)
                else:
                    print("Terminated by user.")
                    return
            else:
                edge_list = load_csv(edge_list_path)
                edge_list1, edge_list2 = get_nearest_neighbors(edge_list,tolerance=1e-6, fraction1=fraction_edge1,fraction2=fraction_edge2)
        else:
            edge_list1 = load_edge_list(edge_list1_path)
            edge_list2 = load_edge_list(edge_list2_path)

        if plot_grid_array:
            # Plot the hexgrid with edges
            plot_hexgrid(hexgrid, edge_list1, edge_list2, show_vertex_numbers=False)
        
        # Merge edge lists
        edge_list_merged = edge_list1 + edge_list2

        # Generate adjacency list and create files
        if run_adjacency_list:
            adjacency_list = generate_adjacency_list(hexgrid, edge_list_merged, fraction_dbonds)
            create_adjacency_list_file(adjacency_list, adjacency_list_path)
            create_output_files(hexgrid, edge_list_merged, output_file)

    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()
