# MDStrainMapper
MDStrainMapper is a comprehensive tool designed to study local deformation in proteins. It provides an effective measure of strain, which captures local deformation relevant to protein function. This repository offers two main functionalities: measuring strain between two configurations and measuring strain over the course of a trajectory.

Table of Contents
=================
* [ Repository Structure](#repository-structure)
* [ Dependencies](#dependencies)
* [ 🚀&nbsp; Getting Started](#🚀&nbsp;-getting-started)
* [ 🧬&nbsp; Usage](#🧬&nbsp;-usage)
  * [ Static Strain](#static-strain)
  * [ Trajectory Processing](#trajectory-processing)
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

## 🧬&nbsp; Usage

### Static Strain 
Measure strain between two configurations (e.g., open and closed) as:

$$
S_i = \frac{1}{N_i} \sum_{\substack{j = 1 \\ j \neq i \\ \Delta_{c}(i, j) < \text{cutoff}}}^{n} \left( \frac{\Delta_{o}(i, j) - \Delta_{c}(i, j)}{\Delta_{c}(i, j)} \right)
$$

##### Step 1: Import Custom Modules 

```python
from src.static_strain import calculate_distance, calculate_strain, plot_strain
```

##### Step 2: Calculate and Plot Strain

```python
from src.trajectory_processing import load_universe

conf1 = load_universe('data/open.pdb')
conf2 = load_universe('data/closed.pdb')

distance_matrix_conf1 = calculate_distance(conf1)
distance_matrix_conf2 = calculate_distance(conf2)

strains = calculate_strain(distance_matrix_conf2, distance_matrix_conf1)
plot_strain(strains, filename='static_strain.png')
```
### Trajectory Processing
Measure strain over the course of a trajectory as:

$$
S_i(t) = \frac{1}{N_i} \sum_{\substack{j = 1 \\ j \neq i \\ }}^{n} \left( \frac{\Delta_{i, j}(t) - \Delta_{i, j}^{\text{initial}}}{\Delta_{i, j}^{\text{initial}}} \right)
for: \Delta_{i, j}^{\text{initial}} < \text{cutoff}
$$

##### Step 1: Import Custom Modules
```python
from src.trajectory_processing import load_universe, process_traj
from src.utils import load_h5_data_to_df, smooth_data, plot_heatmap
```

##### Step 2: Load and Process Trajectory Data
_Note: This step is computationally expensive and requires a GPU._
```python
universe = load_universe('data/1ZNX_topology.gro', '1ZNX_mdtraj.xtc')
process_traj(universe, output='results/local_strain.h5')
```
##### Step 3: Load and Plot Heatmap of Trajectory Data
```python
strain_data_df = load_h5_data_to_df('results/local_strain.h5')
plot_heatmap(strain_data_df, window_size=10, filename='heatmap.png')
```
## Contributing
Contributions are welcome! Please submit a pull request or open an issue to discuss any changes.
