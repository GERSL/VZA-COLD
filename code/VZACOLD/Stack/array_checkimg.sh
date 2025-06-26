#!/bin/bash
#SBATCH -J ImgQA
#SBATCH --partition priority
#SBATCH --ntasks 1
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --constraint='epyc128'
#SBATCH --mem-per-cpu 8G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_checkimg/slurm-%A_%a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_checkimg/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/StackData

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX

export job njobs
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_findimgerror\($job,$njobs\)\;exit\;
