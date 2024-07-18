# MDStrainMapper

MDStrainMapper is a ready-to-use tool to study local deformation in proteins. For a local measure of motion one can examine the displacement per amino acid. However, that displacement can be large in areas that do not undergo any local rearrangements. n alternative approach is to measure strain, which is a local measure of deformation. Strain is low if parts of a protein move together as a rigid body, and high if their relative positions change, which implies stretching, compressing and breaking of bonds, and the formation of new bonds. Thus, strain naturally captures local deformation relevant to function. This tool measures _effective strain_, which is simply the average relative change in distance vectors be- tween neighbors.

This tool can either measure strain between two configuration (eg. open and close) as:

$$
S_i = \frac{1}{N_i} \sum_{\substack{j = 1 \\ j \neq i \\ \Delta_{c}(i, j) < \text{cutoff}}}^{n} \left( \frac{\Delta_{o}(i, j) - \Delta_{c}(i, j)}{\Delta_{c}(i, j)} \right)
$$

For this purpose you will need the pdb structures of the two conformations.

This tool can also measure strain over the course of a trajectory as:

$$
S_i(t) = \frac{1}{N_i} \sum_{\substack{j = 1 \\ j \neq i \\ }}^{n} \left( \frac{\Delta_{i, j}(t) - \Delta_{i, j}^{\text{initial}}}{\Delta_{i, j}^{\text{initial}}} \right)
for: \Delta_{i, j}^{\text{initial}} < \text{cutoff}
$$

For this purpose you will need a trajectory file along with a topology file. 

## How to use

You can either clone the repository directly on a machine with GPU. Note that cuda_cupy12x is required, so insure you have the correct version of CUDA.
You can also use the colab notebook procided or the Docker image. 

## 🚀&nbsp; Getting Started
Step 1: After you clone the repository, setup the dependencies:
```console
pip install -r https://raw.githubusercontent.com/Zahra-Alavi/MDStrainMapper/main/requirements.txt
```

Step 2: Download the Trajectory: for the purpose of this demo, a 50s trajectory is available on Amazon S3 which can be downloaded via:

```console
wget https://mdstrainmapper.s3.us-east-2.amazonaws.com/data/1ZNX_mdtraj.xtc
```

## 📏&nbsp; Static Strain

Step 1: Import Custom Modules:

```python
from src.static_strain import calculate_distance, calculate_strain, plot_strain
```

Step 2: Calculate and plot strain: 

```python
conf1 = load_universe('data/open.pdb')
conf2 = load_universe('data/closed.pdb')

distance_matrix_conf1 = calculate_distance(conf1)
distance_matrix_conf2 = calculate_distance(conf2)

strains = calculate_strain(distance_matrix_conf2, distance_matrix_conf1)
plot_strain(strains, filename='static_strain.png')
```


## 📈&nbsp; Trajectory processing: Strain profile

Step 1: Import Custom Modules:
```python
from src.trajectory_processing import load_universe, process_traj
from src.utils import load_h5_data_to_df, smooth_data, plot_heatmap
```
Step 2: Load and Process Trajectory Data:
```python
# Load the trajectory
universe = load_universe('data/1ZNX_topology.gro', '1ZNX_mdtraj.xtc')
```

_Important_: this step is computationally expensive and requires a GPU. Make sure your CUDA version matches that of the cuda_cupy.

```python
# Process the trajectory and save the output
process_trajectory(universe, output='results/local_strain.h5')
```
Step 3: Load and Plot Heatmap of Trajectory Data
```python
# Load the data
strain_data_df = load_h5_data_to_df('results/local_strain.h5')

# Plot the heatmap
plot_heatmap(df, window_size=10, filename='heatmap.png')
```
