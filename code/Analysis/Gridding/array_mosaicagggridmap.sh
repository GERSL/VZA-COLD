#!/bin/bash
#SBATCH -J MscGrd
#SBATCH -p general
#SBATCH --ntasks 1
#SBATCH --array=1-18
#SBATCH --nodes=1
#SBATCH --constraint epyc128
#SBATCH --mem-per-cpu 32G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_mosaicagggridmap/slurm-%A_%4a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_mosaicagggridmap/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/Analyze/Gridding

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX
export job njobs
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_mosaicagggrid\($job,$njobs\)\;exit\;
