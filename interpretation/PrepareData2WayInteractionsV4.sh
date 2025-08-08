#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --job-name=v4_2way_interactions
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/9.3.0
module load gdal/3.5.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla PrepareData2WayInteractionsV4.R
