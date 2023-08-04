#!/bin/bash -l
#SBATCH -A naiss2023-22-759
#SBATCH -p core
#SBATCH -n 20
#SBATCH -o /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/%x.output
#SBATCH -e /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/%x.output
#SBATCH -t 4-00:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mattias.erhardsson@gu.se

cd /crex/proj/naiss2023-22-759/lmpc-infection-rnaseq/

module load R_packages/4.2.1
Rscript lmpc-infection-rnaseq-analysis-script1-pre-setrank.R --no-restore --no-save
echo R script 1 finished
Rscript lmpc-infection-rnaseq-analysis-script2-setrank.R --no-restore --no-save
echo R script 2 finished, job should now be finished successfully