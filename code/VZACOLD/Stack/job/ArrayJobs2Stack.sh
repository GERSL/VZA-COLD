#!/bin/bash
#SBATCH --job-name JobName
#SBATCH --partition ClusterName
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --array 1-TotalCores
#SBATCH --constraint='epyc128'
#SBATCH --mem=5G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_readstack/slurm-%A_%a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_readstack/error-%A_%4a.err

cd "FolderpathCode"

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX
module load matlab
matlab -nodisplay -singleCompThread -r "loop_createRowdata($job, $njobs, 'TileName', 'FolderpathNTLRaw', 'FolderpathNTLBRDF', 'FolderpathStack'); exit"
