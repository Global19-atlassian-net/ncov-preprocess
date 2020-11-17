#!/bin/sh

#SBATCH --output=logs/cluster_%j.out                 # where to store the output ( %j is the JOBID )
#SBATCH --error=logs/cluster_%j.err                  # where to store error messages

# activate conda environment
source /scicore/home/neher/neher/miniconda3/etc/profile.d/conda.sh
conda activate nextstrain

{exec_job}


