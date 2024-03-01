# Step-by-step MD trajectory analysis tutorial.

## 1. Making Directory and Linking Required Files.
### First, create a new directory named 'analysis' inside your gmx working directory and navigate into it. Then, create symbolic links to the required files from the working directory.

* Create a new directory named 'analysis' and navigate into it.
```
mkdir analysis & cd analysis
```

* Create a symbolic link to the 'md_0_100.tpr' file.
```
ln -s ../md_0_100.tpr
```

* Create a symbolic link to the 'md_0_100.xtc' file.
```
ln -s ../md_0_100.xtc
```
* Create a symbolic link to the 'index.ndx' file.
```
ln -s ../index.ndx
```
## 2. Cleaning Trajectory and Extracting First Frame
### Next, clean the trajectory and extract the first frame using `gmx trjconv`. This tool is used to manipulate and convert trajectory files. The `-pbc`, `-ur`, `-center` and `-n` options are used to remove periodic boundary conditions, make the unit cell compact, center the output and specify the index file, respectively.

* Remove periodic boundary conditions, make the unit cell compact, center the output and specify the index file. Select group for centering: `1` (Protein) and select group for output: `19` (Protein_LIG) when promoted.
```
gmx trjconv -s md_0_100.tpr -f md_0_100.xtc -o md_0_100_noPBC_noWAT.xtc -pbc mol -ur compact -center -n index.ndx
```
* Remove jumps in the trajectory due to periodic boundary conditions. Select the group for centering: `1` (Protein) and select the group for output: `19` (Protein_LIG) when promoted.
```
gmx trjconv -s md_0_100.tpr -f md_0_100_noPBC_noWAT.xtc -o md_0_100_noWAT_noPBC_noJUMP.xtc -pbc nojump -ur compact -center -n index.ndx
```
* Fit the trajectory to a reference structure
```
gmx trjconv -s md_0_100.tpr -f md_0_100_noWAT_noPBC_noJUMP.xtc -o md_0_100_noWAT_noPBC_noJUMP_fit.xtc -ur compact -center -fit rot+trans -n index.ndx
```
* Dump the first frame of the trajectory
```
gmx trjconv -s md_0_100.tpr -f md_0_100_noWAT_noPBC_noJUMP_fit.xtc -o first_frame.pdb -dump 0 -n index.ndx
```
## 3. Essential Structural Analysis: 
### These methods provide insights into the structural changes and stability of the system throughout the simulation. They help in understanding how the structure of the molecule changes over time. This includes Root Mean Square Deviation (RMSD), Root Mean Square Fluctuation (RMSF), Radius of Gyration (ROG/Rg), Hydrogen-bond analysis (HBond), Solvent Accessible Surface Areas (SASA) and Secondary Structure Element (SSE) analysis.

* Create a new directory for RMSD, RMSF, and ROG and navigate into it.
```
mkdir structural_analysis & cd structural_analysis
```
* Calculate the RMSD of the c-alpha atoms from the trajectory
```
gmx rms -s ../md_0_100.tpr -f ../md_0_100_noWAT_noPBC_noJUMP_fit.xtc -o rmsd_ca.xvg -tu ns -fit rot+trans
```
* Calculate the RMSF of the of the c-alpha atoms from the trajectory
```
gmx rmsf -s ../md_0_100.tpr -f ../md_0_100_noWAT_noPBC_noJUMP_fit.xtc -o rmsf_ca.xvg -res -fit
```
* Calculate the ROG of the protein from the trajectory
```
gmx gyrate -s ../md_0_100.tpr -f ../md_0_100_noWAT_noPBC_noJUMP_fit.xtc -o rog_prot.xvg -p
```
* Compute the hydrogen bonds between protein and ligand groups.
```
gmx hbond -s ../md_0_100.tpr -f ../md_0_100_noPBC_noWAT_noJUMP_fit.xtc -n ../index.ndx -num protein-ligand_hbonds.xvg -dan protein-ligand_HBA_HBD.xvg -g hbond.log -dist hbdist.xvg -tu ns -nthreads 70
```
Where, in the `gmx hbond` command:\
`-s md_0_100.tpr`, Specifies the input gmx run topology file, which contains the starting structure of your simulation, the molecular topology, and all the simulation parameters.\
`-f md_0_100_noPBC_noWAT_noJUMP_fit.xtc`, Specifies the input trajectory file, which contains the coordinates of the atoms during the simulation.\
`-n index.ndx`, Specifies the index file, which contains atom groups that you want to analyze.\
`-num protein-ligand_hbonds.xvg`, Specifies the output file containing the number of hydrogen bonds as a function of time.\
`-dan protein-ligand_HBA_HBD.xvg`, Specifies the output file containing information about the average number of donors (HBD) and acceptors (HBA) for each group.\
`-g hbond.log`, Specifies the log file that will contain information about the progress of the calculation.\
`-dist hbdist.xvg`, Specifies the output file containing the distribution of the hydrogen bond lengths.\
`-tu ns`, Specifies the time unit (ns/ps) used in the output files.\
`-nthreads 70`, Specifies the number of threads to be used for the calculation.\

* Comopute the solvent surface area
1. Computes the SASA for Protein.\
```
gmx sasa -s md_0_100.tpr -f md_0_100_noWAT_noPBC_noJUMP_fit.xtc -n index.ndx -o sasa_rec.xvg -or resarea_rec.xvg -surface Protein
```
2. Computes the SASA for Ligand.\
```
gmx sasa -s md_0_100.tpr -f md_0_100_noWAT_noPBC_noJUMP_fit.xtc -n index.ndx -o sasa_lig.xvg -or resarea_lig.xvg -surface TPA
```
3. Computes the SASA for Protein-Ligand Complex.\
```
gmx sasa -s md_0_100.tpr -f md_0_100_noWAT_noPBC_noJUMP_fit.xtc -n index.ndx -o sasa_com.xvg -or resarea_com.xvg -surface Protein_TPA
```
* Computes the number of secondary structure elements for each time frame by calling the dssp program
```
gmx do_dssp -s md_0_100.tpr -f md_0_100_noWAT_noPBC_noJUMP_fit.xtc -tu ns -sss ECBGHITS -o ss.xpm -ssdump ssdump.dat -sc scount.xvg  # Select: Group 1  ( Protein)
```
Where, in `gmx do_dssp` command.\
'E': Parallel Beta-sheet\
'B': Anti-parallel Beta-sheet\
'G': 3-10 helix\
'H': Alpha helix\
'I': Pi (3-14) helix\
'T': Turn\
'S': Bend\
'C': Coil

## 4. Essential Dynamic Analysis:
### These methods provide insights into the dynamic behavior of the system. They help in understanding the movements and interactions of the molecule. This includes PCA (Principal Component Analysis) and Clustering.

### 4.1. PCA: Reduces the dimensionality of the data by identifying the directions (principal components) along which the variation in the data is maximal. It is used to understand the essential dynamics or large-scale motions in the system.
* Create a new directory for PCA and navigate into it
```
mkdir pca_analysis & cd pca_analysis
```
* Activate the appropriate conda environment
```
conda activate mode_task
```
* Perform PCA using 'pca.py' script from MODE-TASK tool.
```
pca.py -p ../first_frame.pdb -t ../md_0_100_noWAT_noPBC_noJUMP_fit.xtc -ag CA -r ../first_frame.pdb -pt svd -out pca 2>&1 | tee pca_log.txt
```
### 4.2. Clustering: Grouping similar structures together from a trajectory helps in understanding the conformational changes during the simulation. Here, we are using the `ttclust` tool for clustering.
* Create a new directory for clustering analysis and navigate into it
```
mkdir clustering_analysis  & cd clustering_analysis  
```
* Activate the conda environment with ttclust installed
```
conda activate ttclust
```
* Run ttclust
```
ttclust -f ../md_0_100_noWAT_noPBC_noJUMP_fit.xtc -t ../first_frame.pdb -l ttclust_clustering.log -st "all" -sa "name CA" -sr "name CA" -m average -aa Y 
```
Where, in the `ttclust` command:\
`-f <traj_filename>`, specifies the input trajectory file (xtc).\
`-t <top_filename>`, specifies the topology file (pdb/gro).\
`-l <out_log_filename>`, specifies the name of out log file.\
`-st "all"`, selects all frames from the trajectory for clustering.\
`-sa "name CA"` and `-sr "name CA"`, select alpha carbons for alignment and RMSD calculation.\
`-m average`, specifies the method of clustering (average linkage).\
`-aa Y`, specifies that all atoms should be used for clustering.


## 5. Binding free energy analysis and its decomposition into the residues of binding site:
### MM/PB(GB)SA analysis is a method to calculate the binding free energy of a system and the decomposition of energy into residues of the active site (within 4/6 Å from the ligand). Here, we are using the `gmx_MMPBSA` tool for this calculation (https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645).

* Create a new directory for MMPBSA analysis and navigate into it
```
mkdir mmpbsa_analysis and cd mmpbsa_analysis
```
* Activate the conda environment with gmx_MMPBSA installed
```
conda activate gmxMMPBSA  
```
* Create input files for MMPBSA calculation and do necessary changes in the parameters accroding to your system.
```
gmx_MMPBSA --create_input pb decomp  
```
* Run gmx_MMPBSA. 
```
mpirun -np 75 gmx_MMPBSA -O -i ../mmpbsa.in -cs ../first_frame.pdb -ci ../index.ndx -cg 1 13 -ct ../md_0_100_noWAT_noPBC_noJUMP_fit.xtc -cp ../../topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -do FINAL_DECOMP_MMPBSA.dat -deo FINAL_DECOMP_MMPBSA.csv -nogui  
```
Where, in gmx_MMPBSA command, Group `1`: Protein and Group `13`: Ligand.
* Analyze the MMPBSA results
```
gmx_MMPBSA_ana  
```
Where, in the `gmx_MMPBSA` command:\
`-O`, Overwrites any existing output files.\
`-i <input_parameter_file>`, specifies the input file for the MMPBSA calculation.\
`-cs <input_structre_filename>`, Specifies the structure file (pdb/gro).\
`-ci <input_index_filename>`, Specifies the index file.\
`-cg 1 13`, Specifies the groups for the calculation.\
`-ct <input_trajectory_filename>`, specifies the trajectory file.\
`-cp <input_topol_filename>`, Specifies the topology file.\
`-o <out_dot_filename>`, Specifies the output dot file for the final results.\
`-eo <out_csv_filename>`,Specifies the output CSV file for the final results.\
`-do <out_dot_filename>`, Specifies the output file for the decomposition.\
`-deo <out_csv_filename>`, Specifies the CSV file for the decomposition.\
`-nogui`, Runs the program without a graphical user interface.
