#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 12:59:29 2024

@author: zalavi
"""

import numpy as np
import pandas as pd

def calculate_distance(universe):
    # Selecting all amino acid residues (adjust the selection if necessary)
    protein = universe.select_atoms('protein')
    n_residues = len(protein.residues)
    distance_matrix = np.zeros((n_residues, n_residues))

    for i, res1 in enumerate(protein.residues):
        for j, res2 in enumerate(protein.residues):
            if i != j:
                # Calculate the distance between the center of masses of the two residues
                distance = np.linalg.norm(res1.atoms.center_of_mass() - res2.atoms.center_of_mass())
                distance_matrix[i, j] = distance
    return distance_matrix

def calculate_strain(matrix_conf1, matrix_conf2, cutoff=15.0, output_csv=None):
    n = matrix_conf1.shape[0]
    strains = np.zeros(n)
    for i in range(n):
        strain_sum = 0
        count = 0
        for j in range(n):
            if i != j and matrix_conf1[i, j] < cutoff:  # Only consider distances within the cutoff
                delta_conf2 = matrix_conf2[i, j]
                delta_conf1 = matrix_conf1[i, j]
                if delta_conf1 != 0:
                    strain_sum += (delta_conf2 - delta_conf1) / delta_conf1
                    count += 1
        if count > 0:
            strains[i] = strain_sum / count  # Average the strain over the number of valid pairs
    if output_csv:     
        # Save strains to CSV
        strain_df = pd.DataFrame(strains, columns=["Strain"])
        strain_df.to_csv(output_csv, index=False)
    
    return strains

