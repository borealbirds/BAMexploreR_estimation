#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=187G
#SBATCH --time=18:00:00
#SBATCH --job-name=v5_cov_importance
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla PrepareDataCovariateImportanceV5.R
