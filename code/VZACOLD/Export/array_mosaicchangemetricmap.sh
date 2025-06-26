#!/bin/bash
#SBATCH -J MscChg
#SBATCH --partition=general
#SBATCH --account=zhz18039
#SBATCH --constraint='epyc128'
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-9
#SBATCH --exclude=cn538
#SBATCH --mem-per-cpu=32G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_mosaicchangemetricmap/%x-%A_%4a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_mosaicchangemetricmap/%x-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/ChangeMap

matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_mosaicchangemetricmap\($SLURM_ARRAY_TASK_ID,$SLURM_ARRAY_TASK_MAX\)\;exit\;
