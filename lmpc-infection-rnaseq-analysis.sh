#!/bin/bash -l
#SBATCH -A naiss2023-22-759
#SBATCH -p node
#SBATCH -N 1
#SBATCH -o /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/%x.output
#SBATCH -e /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/%x.output
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mattias.erhardsson@gu.se

module load R/4.3.1

cd /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/

Rscript 01_lmpc-infection-rnaseq-analysis-environment-setup.R --no-restore --no-save

Rscript 02_lmpc-infection-rnaseq-analysis-pre-setrank.R --no-restore --no-save

Rscript 03_lmpc-infection-rnaseq-analysis-setrank.R --no-restore --no-save