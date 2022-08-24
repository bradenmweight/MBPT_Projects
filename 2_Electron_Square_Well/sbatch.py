#!/software/anaconda3/2020.11/bin/python
#SBATCH -p debug
#SBATCH --job-name=sbatch
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem-per-cpu=2GB
#SBATCH -t 1:00:00
#SBATCH --output=output.out
#SBATCH --error=output.err

import main

main.main()