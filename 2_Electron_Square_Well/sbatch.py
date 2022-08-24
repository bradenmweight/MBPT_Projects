#!/software/anaconda3/2020.11/bin/python
#SBATCH -p debug
#SBATCH --job-name=sbatch
#SBATCH --nodes=1 --ntasks=1
#SBATCH --mem-per-cpu=2GB
#SBATCH -t 1:00:00
#SBATCH --output=output.out
#SBATCH --error=output.err

import os
import sys
sys.path.append(os.popen("pwd").read().replace("\n",""))
import numpy as np

import main

main.main()
#main.get_globals()
#print(main.get_Coulomb_element_single_particle( [1,1,1,1], np.linspace(-5,5,100) ))