# GROMACS Trajectory Analysis Tutorial
This repository contains a step-by-step tutorial for analyzing molecular dynamics (MD) simulations using GROMACS. The tutorial covers the following steps:
*************************************************************************************************************************************************************************************
## 1. Making Directory and Linking Required Files
Create a new directory for analysis and link the required files from the parent directory.

## 2. Cleaning Trajectory and Extracting First Frame
Clean the trajectory and extract the first frame using `gmx trjconv`.

## 3. RMSD, RMSF, ROG/Rg
Calculate the Root Mean Square Deviation (RMSD), Root Mean Square Fluctuation (RMSF), and Radius of Gyration (Rg/ROG) using `gmx rms`, `gmx rmsf`, and `gmx gyrate`, respectively.

## 4. PCA
1. Perform Principal Component Analysis (PCA) using the `pca.py` script from the MODE-TASK tool.
2. Comparing PCA results from two different trajectories

## 5. Clustering
Perform clustering to group similar structures together from a trajectory using the `ttclust` tool.

## 6. MMPBSA
Calculate the binding free energy of a system using the `gmx_MMPBSA` tool.
*************************************************************************************************************************************************************************************

# Detailed instructions and comments for each step are in the tutorial script file.
# Feel free to ask if you have any more questions or need more help. 
# Happy simulating! 😊.



