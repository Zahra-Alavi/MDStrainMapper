#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 12:35:31 2024

@author: zalavi
"""

import h5py
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import MDAnalysis as mda


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


def load_h5_data_to_df(file_path, key='local_strain'):
    with h5py.File(file_path, 'r') as file:
        data = file[key][:]
    return pd.DataFrame(data)

def smooth_data(df, window_size):
    return df.rolling(window=window_size, center=True).mean()


def plot_heatmap(df, window_size=None, filename=None):
    if window_size:
        df = smooth_data(df, window_size)

    plt.figure(figsize=(15, 10))
    sns.heatmap(df, xticklabels=int(df.shape[1] / 10), yticklabels=10, cmap='viridis', cbar_kws={'label': 'Strain'})
    
    # Customize x-axis ticks
    plt.xticks(ticks=np.linspace(0, df.shape[1]-1, 10), labels=np.linspace(0, df.shape[1]-1, 10, dtype=int))
    
    plt.xlabel('Time (Frame Number)')
    plt.ylabel('Residue Number')
    plt.title('Heatmap of Strain vs Time')
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()
 
def plot_strain(strains, filename=None):
    # Calculate the absolute magnitude of strains
    abs_strains = np.abs(strains)
    
    plt.figure(figsize=(10, 5))
    plt.plot(abs_strains, marker='o', linestyle='-', color='b')
    plt.xlabel('Amino Acid Index')
    plt.ylabel('Absolute Strain')
    plt.title('Absolute Strain vs. Amino Acid Number')
    plt.grid(True)
    plt.show()
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()


   
    