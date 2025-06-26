#!/bin/bash
#SBATCH -J IniRad
#SBATCH -p general
##SBATCH --partition priority
##SBATCH --qos=zhz18039epyc
#SBATCH --ntasks 1
#SBATCH --array=1-12
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
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_initialradiance\($job,$njobs\)\;exit\;
