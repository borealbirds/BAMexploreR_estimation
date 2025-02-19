#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --mem=256G
#SBATCH --time=30:00:00
#SBATCH --job-name=v4_cov_importance
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla PrepareDataCovariateImportanceV4.R



#Job ID: 53866707
#Cluster: beluga
#User/Group: mannfred/mannfred
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 25
#CPU Utilized: 04:27:55
#CPU Efficiency: 14.61% of 1-06:33:20 core-walltime
#Job Wall-clock time: 01:13:20
#Memory Utilized: 11.84 GB
#Memory Efficiency: 4.63% of 256.00 GB