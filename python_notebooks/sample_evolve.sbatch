#!/bin/bash
#SBATCH --job-name=make_stream
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nilanjan1@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=hpg2-compute
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=300:00:00
#SBATCH --output=make_stream_%j.out
#SBTACH --error=make_stream_err%j.err
#SBATCH --array=0-200%50

ml conda
source activate galpy1.4
ml gcc gsl python3

python3 sample_evolve_bulge.py --ind $SLURM_ARRAY_TASK_ID
