#!/bin/bash
#SBATCH -p debug
#SBATCH --job-name=sbatch
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem-per-cpu=2GB
#SBATCH -t 1:00:00
#SBATCH --output=main.out
#SBATCH --error=main.err

module load python3
python3 main.py