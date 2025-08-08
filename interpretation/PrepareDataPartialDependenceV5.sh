#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=170G
#SBATCH --time=03:00:00
#SBATCH --job-name=v5_partial_dependence
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla PrepareDataPartialDependenceV5.R


seff 65284974
Job ID: 65284974
Cluster: cedar
User/Group: mannfred/mannfred
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 32
CPU Utilized: 07:15:10
CPU Efficiency: 62.67% of 11:34:24 core-walltime
Job Wall-clock time: 00:21:42
Memory Utilized: 145.66 GB
Memory Efficiency: 85.68% of 170.00 GB (170.00 GB/node)
