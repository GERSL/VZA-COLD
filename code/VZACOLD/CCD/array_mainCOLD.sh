#!/bin/bash
#SBATCH --partition general
##SBATCH --qos=zhz18039epyc
#SBATCH --ntasks 1
#SBATCH -J COLD
#SBATCH --nodes=1
#SBATCH --account=zhz18039
#SBATCH --array=1-3
#SBATCH --constraint epyc128
#SBATCH --mem-per-cpu 16G
#SBATCH --output=/home/til19015/JobOutput/tmp_output_mainCOLD/slurm-%A_%4a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_mainCOLD/error-%A_%4a.err

module purge
module load matlab
cd /home/til19015/GlobalNTLAnalyze/VZACOLD_master/Stable/Map

job=$SLURM_ARRAY_TASK_ID
njobs=$SLURM_ARRAY_TASK_MAX

export job njobs
echo "loop_Convza_irows_v8 $job $njobs"
matlab -nojvm -nodisplay -nosplash -singleCompThread -r loop_Convza_irows\($job,$njobs\)\;exit\;
