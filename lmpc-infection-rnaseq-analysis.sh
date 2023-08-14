#!/bin/bash -l
#SBATCH -A naiss2023-22-759
#SBATCH -p core
#SBATCH -n 16
#SBATCH	-C mem256GB
#SBATCH -o /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/%x.output
#SBATCH -e /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/%x.output
#SBATCH -t 1-10:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mattias.erhardsson@gu.se

module load R/4.3.1

cd /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/

Rscript 01_lmpc-infection-rnaseq-analysis-script1-pre-setrank.R --no-restore --no-save
echo R script 1 finished
Rscript 02_lmpc-infection-rnaseq-analysis-script2-setrank.R --no-restore --no-save
echo R script 2 finished, job should now be finished successfully