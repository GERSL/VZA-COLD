#!/bin/bash
#SBATCH --partition general
#SBATCH --ntasks 1
#SBATCH --array=1-100
#SBATCH --mem-per-cpu 5G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_checkstack/slurm-%A_%a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_checkstack/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/StackData

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX

export job njobs
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_findstackerror\($job,$njobs\)\;exit\;

