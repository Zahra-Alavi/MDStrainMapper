#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:56:35 2024

@author: zalavi
"""

import os
import sys
import argparse
from src.static_strain import calculate_distance, calculate_strain
from src.trajectory_processing import process_traj
from src.utils import load_h5_data_to_df, plot_strain, plot_heatmap, load_universe

def calculate_static_strain(static_conf1_path, static_conf2_path, static_cutoff):
    # Calculate and plot static strain
    conf1 = load_universe(static_conf1_path)
    conf2 = load_universe(static_conf2_path)

    distance_matrix_conf1 = calculate_distance(conf1)
    distance_matrix_conf2 = calculate_distance(conf2)

    strains = calculate_strain(distance_matrix_conf2, distance_matrix_conf1, cutoff=static_cutoff)
    plot_strain(strains, filename='/app/results/static_strain.png')

def calculate_trajectory_strain(topology_path, trajectory_path, trajectory_cutoff, window_size):
    # Load and process trajectory
    universe = load_universe(topology_path, trajectory_path)
    process_traj(universe, output='/app/results/local_strain.h5', cutoff=trajectory_cutoff)

    # Plot strain heatmap
    strain_data_df = load_h5_data_to_df('/app/results/local_strain.h5')
    plot_heatmap(strain_data_df, window_size=window_size, filename='/app/results/heatmap.png')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate protein strain.')
    parser.add_argument('--static', action='store_true', help='Calculate static strain')
    parser.add_argument('--trajectory', action='store_true', help='Calculate trajectory strain')
    parser.add_argument('--static_conf1_path', type=str, help='Path to the first static PDB file')
    parser.add_argument('--static_conf2_path', type=str, help='Path to the second static PDB file')
    parser.add_argument('--topology_path', type=str, help='Path to the topology file')
    parser.add_argument('--trajectory_path', type=str, help='Path to the trajectory file')
    parser.add_argument('--static_cutoff', type=float, default=15.0, help='Cutoff for static strain calculation')
    parser.add_argument('--trajectory_cutoff', type=float, default=15.0, help='Cutoff for trajectory strain calculation')
    parser.add_argument('--window_size', type=int, default=10, help='Window size for smoothing')

    args = parser.parse_args()

    # Ensure results directory exists
    os.makedirs('/app/results', exist_ok=True)

    if args.static:
        if not args.static_conf1_path or not args.static_conf2_path:
            print("Static strain calculation requires --static_conf1_path and --static_conf2_path.")
            sys.exit(1)
        calculate_static_strain(args.static_conf1_path, args.static_conf2_path, args.static_cutoff)
    
    if args.trajectory:
        if not args.topology_path or not args.trajectory_path:
            print("Trajectory strain calculation requires --topology_path and --trajectory_path.")
            sys.exit(1)
        calculate_trajectory_strain(args.topology_path, args.trajectory_path, args.trajectory_cutoff, args.window_size)

    if not args.static and not args.trajectory:
        print("No calculation specified. Use --static and/or --trajectory.")
        sys.exit(1)
