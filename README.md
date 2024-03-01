# **GROMACS Trajectory Analysis Tutorial**
This repository contains a step-by-step tutorial for analyzing molecular dynamics (MD) simulations using GROMACS.

## **MD simulation analysis can be categorized into three main groups:**
### 1. Essential Structural Analysis: 
These methods provide insights into the structural changes and stability of the system throughout the simulation. They help in understanding how the 
structure of the molecule changes over time.
1. **Root Mean Square Deviation (RMSD):** Measures the average distance between the atoms of superimposed proteins.
2. **Root Mean Square Fluctuation (RMSF):** Determines the flexibility of different parts of a protein structure.
3. **Radius of Gyration (ROG/Rg):** This gives an idea about the compactness of the protein structure during the simulation.
4. **Hydrogen-bond Analysis (HBond):** Provides insights into the stability of the protein structure.
5. **Solvent Accessible Surface Areas (SASA):** Calculates the exposure of the protein to the solvent.
6. **Secondary Structure Element (SSE) Analysis:** Provides information about the secondary structure elements of the protein during the simulation.

### 2. Essential Dynamic Analysis:
These methods provide insights into the system's dynamic behaviour. They help to understand the movements and interactions of the molecule.
1. **Principal Component Analysis (PCA).**
2. **Porcupine plot analysis of PCA.**
3. **Clustering analysis.**

### 3. Binding free energy analysis
1. MM-PB(GB)SA
2. Decomposition

## The tutorial covers the following steps:
### 1. Making Directory and Linking Required Files
Create a new directory for analysis and link the required files from the parent directory.

### 2. Cleaning Trajectory and Extracting First Frame
Clean the trajectory and extract the first frame using `gmx trjconv`.

### 3. RMSD, RMSF, ROG/Rg, SASA and SSE
Calculate the Root Mean Square Deviation (RMSD), Root Mean Square Fluctuation (RMSF), and Radius of Gyration (Rg/ROG) using `gmx rms`, `gmx rmsf`, `gmx gyrate`, `gmx sasa`  and `gmx do_dssp` respectively.

### 4. PCA
1. Perform Principal Component Analysis (PCA) using the `pca.py` script from the MODE-TASK tool.
2. Comparing PCA results from two different trajectories

### 5. Clustering
Perform clustering to group similar structures from a trajectory using the `ttclust` tool.

### 6. MMPBSA
Calculate the binding free energy of a system using the `gmx_MMPBSA` tool.
*************************************************************************************************************************************************************************************

# Detailed instructions and comments for each step are in the tutorial script file.
# Feel free to ask if you have any more questions or need more help. 
# Happy simulating! 😊.



