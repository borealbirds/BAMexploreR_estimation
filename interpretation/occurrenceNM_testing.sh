#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --job-name=opticuttest
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/9.3.0
module load gdal/3.9.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla occurrenceNM_testing.R
