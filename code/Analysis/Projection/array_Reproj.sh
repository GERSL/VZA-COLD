#!/bin/bash
#SBATCH -J Reproj
#SBATCH --partition=general
##SBATCH --account=zhz18039
##SBATCH --qos=zhz18039epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-40
#SBATCH --mem-per-cpu=80G
#SBATCH --constraint=epyc128
#SBATCH --output=/home/til19015/JobOutput/tmp_output_projection/%x-%A_%4a.out
#SBATCH --error=/home/til19015/JobOutput/tmp_output_projection/%x-%A_%4a.err

. "/home/til19015/miniconda3_3/etc/profile.d/conda.sh"
conda activate ntlpy310
cd /home/til19015/GlobalNTLAnalyze/Analyze/Projection

python Reproj_Geo2Sin_Tile.py -i $SLURM_ARRAY_TASK_ID -n $SLURM_ARRAY_TASK_MAX
exit
