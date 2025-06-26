#!/bin/bash
#SBATCH -J ExtMap
#SBATCH -p general
#SBATCH --ntasks 1
#SBATCH --array=1-76
#SBATCH --nodes=1
#SBATCH --constraint epyc128
#SBATCH --mem-per-cpu 32G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_extractchangemetricmap/slurm-%A_%4a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_extractchangemetricmap/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/ChangeMap

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX
export job njobs
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_extractdailyabruptprobmap\($job,$njobs\)\;exit\;
