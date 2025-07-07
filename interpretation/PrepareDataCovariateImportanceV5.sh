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



#Job ID: 63606997
#Cluster: cedar
#User/Group: mannfred/mannfred
#State: COMPLETED (exit code 0)
#Nodes: 1
#Cores per node: 32
#CPU Utilized: 05:15:00
#CPU Efficiency: 72.92% of 07:12:00 core-walltime
#Job Wall-clock time: 00:13:30
#Memory Utilized: 160.72 GB
#Memory Efficiency: 85.95% of 187.00 GB (187.00 GB/node)
