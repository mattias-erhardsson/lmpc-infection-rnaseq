#!/bin/bash -l
#SBATCH -A naiss2023-22-759
#SBATCH -p core
#SBATCH -n 16
#SBATCH -C mem256GB
#SBATCH -o ../lmpc-infection-rnaseq/%x.output
#SBATCH -e ../lmpc-infection-rnaseq/%x.output
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mattias.erhardsson@gu.se

module load R/4.3.1

Rscript 01_lmpc-infection-rnaseq-analysis-environment-setup.R --no-restore --no-save

Rscript 02_lmpc-infection-rnaseq-analysis-pre-setrank.R --no-restore --no-save

Rscript 03_lmpc-infection-rnaseq-analysis-setrank.R --no-restore --no-save