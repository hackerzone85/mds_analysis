## Step-by-step MD trajectory analysis tutorial.
1. Making Directory and Linking Required Files
First, create a new directory named 'analysis' inside your gmx working directory and navigate into it. Then, create symbolic links to the required files from the working directory.

mkdir analysis  # Create a new directory named 'analysis'
cd analysis  # Navigate into the 'analysis' directory
ln -s ../md_0_600.tpr .  # Create a symbolic link to the 'md_0_600.tpr' file
ln -s ../md_0_600.xtc .  # Create a symbolic link to the 'md_0_600.xtc' file
ln -s ../index.ndx .  # Create a symbolic link to the 'index.ndx' file

2. Cleaning Trajectory and Extracting First Frame
Next, clean the trajectory and extract the first frame using gmx trjconv. This tool is used to manipulate and convert trajectory files. The -pbc, -ur, -center, and -n options are used to remove periodic boundary conditions, make the unit cell compact, center the output, and specify the index file, respectively.
# Remove periodic boundary conditions, make the unit cell compact, center the output, and specify the index file
gmx trjconv -s md_0_600.tpr -f md_0_600.xtc -o md_0_600_noPBC_noWAT.xtc -pbc mol -ur compact -center -n index.ndx

# Remove jumps in the trajectory due to periodic boundary conditions
gmx trjconv -s md_0_600.tpr -f md_0_600_noPBC_noWAT.xtc -o md_0_600_noPBC_noWAT_noJUMP.xtc -pbc nojump -ur compact -center -n index.ndx

# Fit the trajectory to a reference structure
gmx trjconv -s md_0_600.tpr -f md_0_600_noPBC_noWAT_noJUMP.xtc -o md_0_600_noPBC_noWAT_noJUMP_fit.xtc -ur compact -center -fit rot+trans -n index.ndx

# Dump the first frame of the trajectory
gmx trjconv -s md_0_600.tpr -f md_0_600_noPBC_noWAT_noJUMP_fit.xtc -o first_frame.pdb -dump 0 -n index.ndx

3. RMSD, RMSF, Rg
Calculate the Root Mean Square Deviation (RMSD), Root Mean Square Fluctuation (RMSF), and Radius of Gyration (Rg) using gmx rms, gmx rmsf, and gmx gyrate respectively.
# Calculate the RMSD of the c-alpha atoms from the trajectory
gmx rms -s md_0_600.tpr -f md_0_600_noPBC_noWAT_noJUMP_fit.xtc -o rmsd_ca.xvg -tu ns -fit rot+trans

# Calculate the RMSF of the of the c-alpha atoms from the trajectory
gmx rmsf -s md_0_600.tpr -f md_0_600_noPBC_noWAT_noJUMP_fit.xtc -o rmsf_ca.xvg -res -fit

# Calculate the radius of gyration of the protein from the trajectory
gmx gyrate -s md_0_600.tpr -f md_0_600_noPBC_noWAT_noJUMP_fit.xtc -o rog_prot.xvg -p

4. PCA
Finally, perform Principal Component Analysis (PCA) using the pca.py script. Make sure to activate the appropriate conda environment before running this step.
# Activate the appropriate conda environment
conda activate mode_task

# Create a new directory for PCA and navigate into it
mkdir pca
cd pca

# Perform PCA
pca.py -p ../first_frame.pdb -t ../md_0_600_noPBC_noWAT_noJUMP_fit.xtc -ag CA -r ../first_frame.pdb -pt svd -out pca

5. Clustering
Clustering is a method used to group similar structures together from a trajectory. It helps in understanding the conformational changes during the simulation. Here, we are using the `ttclust` tool for clustering.

conda activate ttclust  # Activate the conda environment with ttclust installed
mkdir clustering  # Create a new directory for clustering analysis
cd clustering  # Navigate into the 'clustering' directory
ttclust -f ../md_0_600_noPBC_noWAT_noJUMP_fit.xtc -t ../first_frame.pdb -l ttclust_clustering.log -st "all" -sa "name CA" -sr "name CA" -m average -aa Y  # Run ttclust

In the `ttclust` command:
- `-f` specifies the input trajectory file.
- `-t` specifies the topology file.
- `-l` specifies the log file.
- `-st "all"` selects all frames from the trajectory for clustering.
- `-sa "name CA"` and `-sr "name CA"` select alpha carbons for alignment and RMSD calculation.
- `-m average` specifies the method of clustering (average linkage).
- `-aa Y` specifies that all atoms should be used for clustering.

6. MMPBSA
MMPBSA (Molecular Mechanics Poisson-Boltzmann Surface Area) is a method to calculate the binding free energy of a system. Here, we are using the `gmx_MMPBSA` tool for this calculation.

conda activate gmxMMPBSA  # Activate the conda environment with gmx_MMPBSA installed
gmx_MMPBSA --create_input pb decomp  # Create input files for MMPBSA calculation
mkdir mmpbsa_analysis  # Create a new directory for MMPBSA analysis
cd mmpbsa_analysis  # Navigate into the 'mmpbsa_analysis' directory
mpirun -np 75 gmx_MMPBSA -O -i ../mmpbsa.in -cs ../first_frame.pdb -ci ../index.ndx -cg 1 13 -ct ../md_0_600_noPBC_noWAT_noJUMP_fit.xtc -cp ../../topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -do FINAL_DECOMP_MMPBSA.dat -deo FINAL_DECOMP_MMPBSA.csv -nogui  # Run gmx_MMPBSA
gmx_MMPBSA_ana  # Analyze the MMPBSA results

In the `gmx_MMPBSA` command:
- `-O` overwrites any existing output files.
- `-i` specifies the input file for the MMPBSA calculation.
- `-cs` specifies the structure file.
- `-ci` specifies the index file.
- `-cg 1 13` specifies the groups for the calculation.
- `-ct` specifies the trajectory file.
- `-cp` specifies the topology file.
- `-o` specifies the output file for the final results.
- `-eo` specifies the CSV file for the final results.
- `-do` specifies the output file for the decomposition.
- `-deo` specifies the CSV file for the decomposition.
- `-nogui` runs the program without a graphical user interface.
