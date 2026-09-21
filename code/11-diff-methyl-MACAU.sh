#!/bin/bash
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=macau_%j.out
#SBATCH --error=macau_%j.err
#SBATCH --chdir=/gscratch/srlab/kdurkin1/ceasmallr/code

apptainer exec --bind /gscratch:/gscratch \
 /gscratch/srlab/kdurkin1/srlab-R4.4-bioinformatics-container-703094b.sif \
  bash -c 'source /srlab/programs/miniforge3-24.7.1-0/etc/profile.d/conda.sh && conda activate /gscratch/srlab/kdurkin1/.conda/envs/macau && Rscript -e "rmarkdown::render(\"11-diff-methyl-MACAU.Rmd\")"'
