#!/bin/bash
#SBATCH -J ChgMtr
#SBATCH -p general
#SBATCH --ntasks 1
#SBATCH --array=1-500
#SBATCH --nodes=1
#SBATCH --constraint epyc128
#SBATCH --mem-per-cpu 5G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_changemetricmap/slurm-%A_%4a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_changemetricmap/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/ChangeMap

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX

export job njobs
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_changemetricmap_irows_ConCOLD_day_Lebanon\($job,$njobs\)\;exit\;
