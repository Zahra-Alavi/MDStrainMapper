#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:23:54 2024

@author: zalavi
"""

import h5py
import pandas as pd

def load_h5_data_to_df(file_path, key='changes_sum'):
    with h5py.File(file_path, 'r') as file:
        data = file[key][:]
    return pd.DataFrame(data)

def smooth_data(df, window_size):
    return df.rolling(window=window_size, center=True).mean()
