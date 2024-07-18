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

def load_h5_data_to_df(file_path, key='local_strain'):
    with h5py.File(file_path, 'r') as file:
        data = file[key][:]
    return pd.DataFrame(data)

def smooth_data(df, window_size):
    return df.rolling(window=window_size, center=True).mean()


def plot_heatmap(df, window_size=None, filename=None):
    
    if window_size:
        df = smooth_data(df, window_size)
        
    # Convert x-axis values to ns (steps multiplied by 0.0125)
    time_ns = df.columns * 0.0125
    df.columns = time_ns
    
    plt.figure(figsize=(15, 10))
    sns.heatmap(df, xticklabels=50, yticklabels=10, cmap='viridis', cbar_kws={'label': 'Strain'})
    plt.xlabel('Time (ns)')
    plt.ylabel('Residue Number')
    plt.title('Heatmap of Strain vs Time')
    
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()
    
    
    