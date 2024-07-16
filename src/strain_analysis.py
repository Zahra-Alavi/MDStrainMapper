#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:21:06 2024

@author: zalavi
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils import load_h5_data_to_df, smooth_data

def calculate_distance_matrix(universe, selection='protein'):
    protein = universe.select_atoms(selection)
    n_residues = len(protein.residues)
    distance_matrix = np.zeros((n_residues, n_residues))
    for i, res1 in enumerate(protein.residues):
        for j, res2 in enumerate(protein.residues):
            if i != j:
                distance = np.linalg.norm(res1.atoms.center_of_mass() - res2.atoms.center_of_mass())
                distance_matrix[i, j] = distance
    return distance_matrix

def calculate_strain(matrix_c, matrix_o, cutoff=15.0):
    n = matrix_c.shape[0]
    strains = np.zeros(n)
    for i in range(n):
        strain_sum = 0
        count = 0
        for j in range(n):
            if i != j and matrix_c[i, j] < cutoff:
                delta = matrix_o[i, j] - matrix_c[i, j]
                if matrix_c[i, j] != 0:
                    strain_sum += delta / matrix_c[i, j]
                    count += 1
        if count > 0:
            strains[i] = strain_sum / count
    return strains

def plot_strain(strains, title='Strain vs. Amino Acid Index', filename=None):
    plt.figure(figsize=(10, 5))
    plt.plot(np.abs(strains), marker='o', linestyle='-', color='b')
    plt.xlabel('Amino Acid Index')
    plt.ylabel('Absolute Strain')
    plt.title(title)
    plt.grid(True)
    if filename:
        plt.savefig(filename)
    else:
        plt.show()
    plt.close()

def plot_heatmap(df, title, window_size=None, filename=None):
    if window_size:
        df = smooth_data(df, window_size)
    sns.heatmap(df, cmap='viridis', cbar_kws={'label': 'Strain'})
    plt.xlabel('Time (ns)')
    plt.ylabel('Residue Number')
    plt.title(title)
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()
