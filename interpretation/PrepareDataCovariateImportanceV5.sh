#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=187G
#SBATCH --time=24:00:00
#SBATCH --job-name=v5_cov_importance
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla PrepareDataCovariateImportanceV5.R


Job ID: 64926931
Cluster: cedar
User/Group: mannfred/mannfred
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 32
CPU Utilized: 07:32:53
CPU Efficiency: 70.35% of 10:43:44 core-walltime
Job Wall-clock time: 00:20:07
Memory Utilized: 154.21 GB
Memory Efficiency: 82.47% of 187.00 GB (187.00 GB/node)