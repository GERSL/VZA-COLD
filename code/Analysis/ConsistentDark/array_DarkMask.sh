#!/bin/bash
#SBATCH -J DrkMsk
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=5G
#SBATCH --constraint='epyc128'
#SBATCH --output=log/%x-%A_%4a.out
#SBATCH --error=log/%x-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/Analyze/DarkPixel

matlab -nosplash -singleCompThread -r loop_DarkMask\($SLURM_ARRAY_TASK_ID,$SLURM_ARRAY_TASK_MAX\)\;exit\;
