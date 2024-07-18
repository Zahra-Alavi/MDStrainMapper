#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:23:36 2024

@author: zalavi
"""

import MDAnalysis as mda
import cupy as cp
import h5py

def load_universe(trajectory_file, topology_file=None):
    """
    Load the MDAnalysis Universe from a topology and trajectory file.
    :param trajectory_file: Path to the trajectory file.
    :param topology_file: Path to the topology file. Optional if trajectory file contains topology information.
    :return: MDAnalysis Universe object.
    """
    try:
        if topology_file:
            universe = mda.Universe(topology_file, trajectory_file)
        else:
            universe = mda.Universe(trajectory_file)
        print("Universe loaded successfully")
        print("Number of atoms:", len(universe.atoms))
        print("Total number of residues:", len(universe.residues))
        return universe
    except Exception as e:
        print("Failed to load universe. Error:", str(e))
        return None  # or raise an exception depending on how you want to handle errors

def process_trajectory(universe, output_file='changes_sum.h5'):
    protein = universe.select_atoms('protein')
    num_residues = len(protein.residues)
    all_centers_of_mass = cp.array([res.atoms.center_of_mass() for res in protein.residues])
    initial_distances = cp.linalg.norm(cp.expand_dims(all_centers_of_mass, axis=0) - cp.expand_dims(all_centers_of_mass, axis=1), axis=-1)
    cp.fill_diagonal(initial_distances, 0)

    changes_sum = cp.zeros((num_residues, len(universe.trajectory)))

    for ts in universe.trajectory:
        current_centers_of_mass = cp.array([res.atoms.center_of_mass() for res in protein.residues])
        current_distances = cp.linalg.norm(cp.expand_dims(current_centers_of_mass, axis=0) - cp.expand_dims(current_centers_of_mass, axis=1), axis=-1)
        cp.fill_diagonal(current_distances, 0)
        strain = (current_distances - initial_distances) / initial_distances
        changes_sum[:, ts.frame] = cp.nansum(strain, axis=1) / cp.sum(strain != 0, axis=1)

    changes_sum_np = cp.asnumpy(changes_sum)
    with h5py.File(output_file, 'w') as f:
        f.create_dataset('changes_sum', data=changes_sum_np)
