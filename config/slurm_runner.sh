#!/bin/zsh

#SBATCH --job-name=euclidlike_sim
#SBATCH --partition=RM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/path/to/slurm-%A-%a.out
#SBATCH --error=/path/to/slurm-%A-%a.err
#SBATCH --array=0-576%200
###SBATCH --array=0-0
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --mem=8G

source activate [env_name]

# Check we are in the right environment
echo $CONDA_PREFIX

export OMP_NUM_THREADS=1

file_list='/path/to/run_list.txt'

# Check we use the right run list
echo ${file_list}

# The following line read the run_list.txt and get the pointing and CCD_ID
i=$((SLURM_ARRAY_TASK_ID + 1))

line=$(sed "${i}q;d" "$file_list")

# Extract both columns using awk
first_column=$(echo "$line" | awk '{print $1}')
second_column=$(echo "$line" | awk '{print $2}')

# Run imsim
galsim was.yaml input.obseq_data.visit=${first_column} image.CCD=${second_column}
