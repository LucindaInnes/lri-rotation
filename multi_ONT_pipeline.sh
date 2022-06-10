#!/bin/bash -l

#SBATCH --partition=panda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=mapping
#SBATCH --time=24:00:00
#SBATCH --mem=32G

source ~/.bashrc

cd $1

for DIR in */; do sbatch $2 $DIR $3; done
