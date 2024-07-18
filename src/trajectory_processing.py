#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:35:10 2024

@author: zalavi
"""

import MDAnalysis as mda
import cupy as cp
import numpy as np
import h5py

def load_universe(trajectory, topology = None):
    """
    Load the MDAnalysis Universe from a topology and trajectory file.
    trajectory: Path to the trajectory file.
    topology: Path to the topology file. Optional if trajectory file contains topology information.
    return: MDAnalysis Universe object.
    """
    if topology:
            universe = mda.Universe(topology, trajectory)
    else:
            universe = mda.Universe(trajectory)
    print("Universe loaded successfully")
    print("Number of atoms:", len(universe.atoms))
    print("Total number of residues:", len((universe.select_atoms('protein')).residues))
    return universe

    

def process_traj(universe, output = 'strain.h5'):
    protein = universe.select_atoms('protein')
    num_residues = len(protein.residues)
    # Initialize GPU array for centers of mass for all residues
    all_centers_of_mass = cp.array([res.atoms.center_of_mass() for res in protein.residues])

    # Calculate initial distances using efficient GPU operations
    initial_distances = cp.linalg.norm(cp.expand_dims(all_centers_of_mass, axis=0) - cp.expand_dims(all_centers_of_mass, axis=1), axis=-1)
    cp.fill_diagonal(initial_distances, 0)

    # Apply cutoff for distances
    cut_off = 1500.0
    valid_distances = initial_distances < cut_off

    print(f"Total number of frames in the trajectory: {len(universe.trajectory)}")

    local_strain = cp.zeros((num_residues, len(universe.trajectory)))


    for ts in universe.trajectory:
        if ts.time % 1000 == 0:
            print(f"Processing frame {ts.frame} at time {ts.time} ps")
    
        # Compute centers of mass for the current frame
        current_centers_of_mass = cp.array([res.atoms.center_of_mass() for res in protein.residues])
        current_distances = cp.linalg.norm(cp.expand_dims(current_centers_of_mass, axis=0) - cp.expand_dims(current_centers_of_mass, axis=1), axis=-1)
        cp.fill_diagonal(current_distances, 0)

        # Avoid zero distances in current distances as well
        current_distances = cp.where(current_distances == 0, cp.nan, current_distances)

        # Compute strain normalized by initial distances where distances are less than the cutoff
        initial_distances_nonzero = cp.where(valid_distances, initial_distances, 1)  # Avoid division by zero
        strain = cp.where(valid_distances, cp.abs(current_distances - initial_distances) / initial_distances_nonzero, 0)
        local_strain[:, ts.frame] = cp.nansum(strain, axis=1) / cp.sum(valid_distances, axis=1)

        print("Finished processing")
        print("Saving data using h5py")

        # Save the data
        local_strain_np = cp.asnumpy(local_strain)
        with h5py.File('local_strain3.h5', 'w') as f:
            f.create_dataset('local_strain', data=local_strain_np)