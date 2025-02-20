#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --job-name=v5_cov_importance
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla PrepareDataCovariateImportanceV5.R


# seff 53897245
# Job ID: 53897245
# Cluster: beluga
# User/Group: mannfred/mannfred
# State: FAILED (exit code 1) (Later finished by running the second half of the script with usage outlined below)
# Nodes: 1
# Cores per node: 25
# CPU Utilized: 05:55:29
# CPU Efficiency: 20.08% of 1-05:30:25 core-walltime
# Job Wall-clock time: 01:10:49
# Memory Utilized: 38.83 GB
# Memory Efficiency: 60.68% of 64.00 GB

# Job ID: 53928834 (second half after fixing an error)
# Cluster: beluga
# User/Group: mannfred/mannfred
# State: COMPLETED (exit code 0)
# Cores: 1
# CPU Utilized: 00:34:17
# CPU Efficiency: 98.09% of 00:34:57 core-walltime
# Job Wall-clock time: 00:34:57
# Memory Utilized: 429.76 MB
# Memory Efficiency: 1.31% of 32.00 GB