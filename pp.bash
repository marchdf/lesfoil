#!/usr/bin/env bash

#SBATCH --job-name=pp
#SBATCH --account=exawind
#SBATCH --nodes=2
#SBATCH --time=04:00:00
#SBATCH -o %x.o%j

export MPI_IGNORE_PBS=on
source /nopt/nrel/ecom/exawind/exawind/scripts/exawind-env-python.sh

mpiexec -np 216 python pp.py -m /scratch/jmelvin/newNalu/baseNalu/forMarc/test/nalu-wind/AAIRFOIL/results/CO2_lesfoil-fine-altSetup-ams-rst2.e
