#!/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=emma.corre@inrae.fr
#SBATCH --mail-type=ALL
#SBATCH --job-name=minimap

source /etc/profile.d/modules.sh
module load minimap2/2.24

genome=$1
longread=$2
repo=$(pwd)

minimap2 -ax asm20 ${genome} ${longreads} --secondary=no -o mapping.sam
