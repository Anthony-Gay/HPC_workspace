#!/bin/bash
#
#SBATCH --ntasks=4
#SBATCH --job-name=foam_test
#SBATCH --output=output.%j.foam_test
#SBATCH --begin=now
#SBATCH --time=00:05:00

####

module add OpenFOAM/7
./Allrun