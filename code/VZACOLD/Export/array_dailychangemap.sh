#!/bin/bash
#SBATCH --partition priority
#SBATCH --ntasks 1
#SBATCH --array=1-120
#SBATCH --constraint epyc128
#SBATCH --mem-per-cpu 6G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_dailychangemap/slurm-%A_%4a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_dailychangemap/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/ChangeMap

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX
export job njobs
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_dailychangemap\($job,$njobs\)\;exit\;
