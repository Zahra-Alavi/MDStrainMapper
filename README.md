# MDStrainMapper
MDStrainMapper is a comprehensive tool designed to study local deformation in proteins. It provides an effective measure of strain, which captures local deformation relevant to protein function [John McBride & Tsvi Tlusty (2024)](https://doi.org/10.1088/1742-5468/ad1be7). This repository offers two main functionalities: measuring strain between two configurations and measuring strain over the course of a trajectory.

Table of Contents
=================
* [ Repository Structure](#repository-structure)
* [ Dependencies](#dependencies)
* [ 🚀&nbsp; Getting Started](#install)
* [ 🧬&nbsp; Usage](#usage)
  * [ Static Strain](#static-strain)
  * [ Trajectory Processing](#trajectory-processing)
* [ 🐳&nbsp; Using the Docker Image with Your Own Data](#docker)
* [ Contributing](#contributing)


## Repository Structure
- `data/`: Contains example PDB and topology files.
- `src/`: Contains the Python modules for strain calculation and trajectory processing.
- `requirements.txt`: Lists the dependencies required to run the package.
- `README.md`: Provides an overview and instructions for using the package.

## Dependencies
MDStrainMapper requires the following Python packages:

- cupy-cuda12x==13.2.0
- h5py==3.11.0
- matplotlib==3.8.4
- MDAnalysis==2.7.0
- numpy==1.26.4
- pandas==2.2.2
- seaborn==0.13.2

Ensure you have the correct version of CUDA to use CuPy effectively.

<a name="install"></a>
## 🚀&nbsp; Getting Started

##### Step 1: Clone the Repository

```sh
git clone https://github.com/Zahra-Alavi/MDStrainMapper.git
cd MDStrainMapper
```
##### Step 2: Install Dependencies
```sh
pip install -r requirements.txt
```
##### Step 3: Download the Trajectory File

For demonstration purposes, download a sample trajectory file: 
```sh 
wget https://mdstrainmapper.s3.us-east-2.amazonaws.com/data/1ZNX_mdtraj.xtc
```
<a name="usage"></a>
## 🧬&nbsp; Usage

### Static Strain 
Measures strain of each amino acid, between two configurations (e.g., open and closed) as:

$$
S_i = \frac{1}{N_i} \sum_{\substack{j = 1 \\ j \neq i \\ \Delta_{c}(i, j) < \text{cutoff}}}^{n} \left( \frac{\Delta_{o}(i, j) - \Delta_{c}(i, j)}{\Delta_{c}(i, j)} \right)
$$

where $S_i$ is the net strain (relative distance change) on residue i, $\Delta_{c}(i, j)$ is the distance between residues i and j in the close conformation and $\Delta_{o}(i, j)$ is the distance in the open conformation. The sum is over all the neighbors as defined by the cutoff distance. $N_i$ is the number of neighbor residues (residues within the cutoff distance.)

##### Step 1: Import Custom Modules 

```python
from src.utils import load_universe, plot_strain
from src.static_strain import calculate_distance, calculate_strain
```

##### Step 2: Calculate and Plot Strain

```python

conf1 = load_universe('data/open.pdb')
conf2 = load_universe('data/closed.pdb')

distance_matrix_conf1 = calculate_distance(conf1)
distance_matrix_conf2 = calculate_distance(conf2)

strains = calculate_strain(distance_matrix_conf2, distance_matrix_conf1, output_csv='static_strain.csv')
plot_strain(strains, filename='static_strain.png')
```
### Trajectory Processing
Measure strain over the course of a trajectory as:

$$
S_i(t) = \frac{1}{N_i} \sum_{\substack{j = 1 \\ j \neq i \\ }}^{n} \left( \frac{\Delta_{i, j}(t) - \Delta_{i, j}^{\text{initial}}}{\Delta_{i, j}^{\text{initial}}} \right)
for: \Delta_{i, j}^{\text{initial}} < \text{cutoff}
$$

where $S_i(t)$ is the net strain (relative distance change) on residue i at time t, $\Delta_{i, j}^{\text{initial}}$ is the distance between residues i and j at the beginning of the simulation ($t=0$) and $\Delta_{i, j}(t)$ is the distance between residues i and j at time t. The sum is over all the neighbors as defined by the cutoff distance. $N_i$ is the number of neighbor residues (residues within the cutoff distance.)


##### Step 1: Import Custom Modules
```python
from src.utils import load_h5_data_to_df, smooth_data, plot_heatmap, load_universe
from src.trajectory_processing import process_traj
```

##### Step 2: Load and Process Trajectory Data
_Note: This step is computationally expensive and requires a GPU._
```python
universe = load_universe('1ZNX_mdtraj.xtc', 'data/1ZNX_topology.gro')
process_traj(universe, output='results/local_strain.h5')
```
##### Step 3: Load and Plot Heatmap of Trajectory Data
```python
strain_data_df = load_h5_data_to_df('results/local_strain.h5')
plot_heatmap(strain_data_df, window_size=10, filename='heatmap.png')
```
<a name="docker"></a>
## 🐳&nbsp; Using the Docker Image with Your Own Data

To use the MDStrainMapper Docker image with your own static PDB files, trajectory data, and custom parameters, follow these steps:


### Case 1: Calculate Static Strain Only

1. **Prepare Your Data**:
   - Ensure you have your static PDB files ready.
   - Place them in a directory on your host machine.

2. **Pull the Docker Image**:
   - Run the following command to pull the Docker image from GitHub Container Registry:

   ```sh
   docker pull ghcr.io/zahra-alavi/mdstrainmapper:latest
   ```

3. **Run the Docker Container**:
   - Use the following command to run the Docker container, replacing the paths to your own data files and setting your custom parameters:

   ```sh
   docker run --gpus all -v /path/to/your/data:/app/data -v /path/to/your/results:/app/results ghcr.io/zahra-alavi/mdstrainmapper:latest \
    --static --static_conf1_path /app/data/conf1.pdb --static_conf2_path /app/data/conf2.pdb --static_cutoff 15.0
   ```
### Case 2: Calculate Trajectory Strain Only (GPU needed)

1. **Prepare Your Data**:
   - Ensure you have your topology and trajectory files ready.
   - Place them in a directory on your host machine.

2. **Pull the Docker Image**:
   - Run the following command to pull the Docker image from GitHub Container Registry:

   ```sh
   docker pull ghcr.io/zahra-alavi/mdstrainmapper:latest
   ```
3. **Run the Docker Container**:
   - Use the following command to run the Docker container, replacing the paths to your own data files and setting your custom parameters:
  
   ```sh
   docker run --gpus all -v /path/to/your/data:/app/data -v /path/to/your/results:/app/results ghcr.io/zahra-alavi/mdstrainmapper:latest \
    --trajectory --topology_path /app/data/topology.gro --trajectory_path /app/data/trajectory.xtc --trajectory_cutoff 15.0 --window_size 10
   ```
### Case 3: Calculate Both Static and Trajectory Strain (GPU needed)

1. **Prepare Your Data**:
   - Ensure you have your static PDB files, topology, and trajectory files ready.
   - Place them in a directory on your host machine.

2. **Pull the Docker Image**:
   - Run the following command to pull the Docker image from GitHub Container Registry:

   ```sh
   docker pull ghcr.io/zahra-alavi/mdstrainmapper:latest
   ```
3. **Run the Docker Container**:
   - Use the following command to run the Docker container, replacing the paths to your own data files and setting your custom parameters:
   ```sh
   docker run --gpus all -v /path/to/your/data:/app/data -v /path/to/your/results:/app/results ghcr.io/zahra-alavi/mdstrainmapper:latest \
    --static --static_conf1_path /app/data/conf1.pdb --static_conf2_path /app/data/conf2.pdb --static_cutoff 15.0 \
    --trajectory --topology_path /app/data/topology.gro --trajectory_path /app/data/trajectory.xtc --trajectory_cutoff 15.0 --window_size 10
   ```
## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.
