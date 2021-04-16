#!/bin/bash -l
#SBATCH --account=harmslab      ### change this to your actual account for charging
#SBATCH --job-name=extract      ### job name
#SBATCH --output=hostname.out   ### file in which to store job stdout
#SBATCH --error=hostname.err    ### file in which to store job stderr
#SBATCH --partition=short  
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1 

# Recommended workflow from gromacs people
# 1. First make your molecules whole if you want them whole.
# 2. Cluster your molecules/particles if you want them clustered.
# 3. Extract the first frame from the trajectory as reference for removing jumps
#    if you want to remove jumps.
# 4. Remove jumps if you want to have them removed using the first frame
# 5. Center your system using some criterion. Doing so shifts the system, so
#    don't use trjconv -pbc nojump after this step.
# 6. Put everything in some box.
# 7. Fit if desired and don't use any PBC related option afterwards.

module load intel/19
module load gromacs/2019.4

# Parse command line, checking for errors.
USAGE="extract-traj.sh base_name starting_gro dest_directory [split_interval_ps]"

base_name=${1}
starting_gro=${2}
dest_dir=${3}
if [ ! "${base_name}" ] || [ ! "${starting_gro}" ] || [ ! "${dest_dir}" ]; then
    echo $USAGE
    exit
fi
split_interval=${4}

xtc_file="${base_name}.xtc"
tpr_file="${base_name}.tpr"
out_root="traj"

# Nuke output directory, remake, and copy in starting_gro file
rm -rf ${dest_dir}
mkdir ${dest_dir}

# Copy the gro file to the output folder
cp ${starting_gro} ${dest_dir}/${out_root}.gro
cp ${tpr_file} ${dest_dir}/${out_root}.tpr

# Use trjconv to write out all protein and solvent atoms
mkdir extract-tmp
cd extract-tmp

# Make whole 
echo "Make whole."
echo "0" > junk.command
cat junk.command | gmx trjconv -f ../${xtc_file} -pbc whole -o ${out_root}.xtc -s ../${tpr_file}

# Now remove jumps
echo "Remove jumps."
echo "0" > junk.command
cat junk.command | gmx trjconv -f ${out_root}.xtc -pbc nojump -o ${out_root}.xtc -s ../${tpr_file}

# Now center on the protein
echo "Center on protein."
echo "1" > junk.command
echo "0" >> junk.command
cat junk.command | gmx trjconv -f ${out_root}.xtc -center -pbc mol -o ${out_root}.xtc -s ../${tpr_file}

# Now align
echo "3" > junk.command
echo "0" >> junk.command
cat junk.command | gmx trjconv -f ${out_root}.xtc -s ../${tpr_file} -fit rot+trans -o ${out_root}.xtc

# Split the final trajectory if desired, then move to destination directory
if [ ${split_interval} ]; then
    echo "Split file."
    echo "0" > junk.command
    cat junk.command | gmx trjconv -f ${out_root}.xtc -split ${split_interval} -s ../${tpr_file} -o ${out_root}-split.xtc 

    cp *split* ../${dest_dir}
else
    cp ${out_root}.xtc ../${dest_dir}
fi

# Delete extraction directory
cd ..
rm -rf extract-tmp

