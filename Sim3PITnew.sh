#!/bin/bash

#SBATCH -J SimMultiBaseIndAccess
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -p longjobs
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH -t 14-00:00:00
#SBATCH --mail-user aramir21@eafit.edu.co
#SBATCH --mail-type ALL

module load python/3.6.0_miniconda-4.3.11_gcc-11.2.0
source activate ARH

Rscript SimulationV2C.R
