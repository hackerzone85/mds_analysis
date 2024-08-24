#!/bin/bash

# Define the systems to process
systems=("6_MPro_Wild_PLZ" "7_MPro_Mute_PLZ" "8_MPro_Wild_BTZ" "9_MPro_Mute_BTZ" "10_MPro_Wild_OGT" "11_MPro_Mute_OGT") 

# Load the gmx command
source /home/ddalab/GMX2020.3/bin/GMXRC

# ANSI escape code for yellow
YELLOW='\033[1;33m'
NOCOLOR='\033[0m'

# Loop over each system
for sys in "${systems[@]}"; do
    echo -e "${YELLOW}Processing system: $sys${NOCOLOR}" | boxes -d info
    cd "$sys" || { echo -e "${YELLOW}Failed to enter directory $sys${NOCOLOR}" | boxes -d WARNING -p a2v1; exit 1; }

    # Step 1: Check if md production job was previously run
    if [ -f "md_0_100.cpt" ]; then
        echo -e "${YELLOW}Step 1: md_0_100.cpt exists. Resuming the mdrun job...${NOCOLOR}" | boxes -d info -p a2v1
        gmx mdrun -ntmpi 12 -ntomp 6 -pin on -v -nb gpu -pme cpu -tunepme yes -pmefft cpu -bonded gpu -update cpu -deffnm md_0_100 -cpi md_0_100.cpt -s md_0_100.tpr -append
    else
        # Step 2: Check for the existence of md_0_100.tpr file
        if [ -f "md_0_100.tpr" ]; then
            echo -e "${YELLOW}Step 2: md_0_100.tpr exists. Starting the mdrun job...${NOCOLOR}" | boxes -d info -p a2v1
            gmx mdrun -ntmpi 12 -ntomp 6 -pin on -v -nb gpu -pme cpu -tunepme yes -pmefft cpu -bonded gpu -update cpu -deffnm md_0_100
        else
            # Step 3: Generate the md_0_100.tpr file since it does not exist
            echo -e "${YELLOW}Step 3: md_0_100.tpr does not exist. Generating it now...${NOCOLOR}" | boxes -d info -p a2v1
            gmx grompp -v -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -po mdout.mdp -pp processed.top -o md_0_100.tpr
            
            # Step 4: Run the md production job after generating the tpr file
            echo -e "${YELLOW}Step 4: Running the mdrun job after generating md_0_100.tpr...${NOCOLOR}" | boxes -d info -p a2v1
            gmx mdrun -ntmpi 12 -ntomp 6 -pin on -v -nb gpu -pme cpu -tunepme yes -pmefft cpu -bonded gpu -update cpu -deffnm md_0_100
        fi
    fi

    # Return to the parent directory
    cd .. || { echo -e "${YELLOW}Failed to return to the parent directory${NOCOLOR}" | boxes -d WARNING -p a2v1; exit 1; }
done
echo -e "${YELLOW}All systems processed successfully!${NOCOLOR}" | boxes -d info

