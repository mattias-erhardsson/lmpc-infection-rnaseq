#!/bin/bash -l
#SBATCH -A naiss2023-22-759
#SBATCH -p core
#SBATCH -n 20
#SBATCH -o /crex/proj/naiss2023-22-759/%x.output
#SBATCH -e /crex/proj/naiss2023-22-759/%x.output
#SBATCH -t 4-00:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=mattias.erhardsson@gu.se

cd /crex/proj/naiss2023-22-759/

module load R_packages/4.2.1
module load RStudio/2022.07.1-554
Rscript lmpc-infection-rnaseq-analysis.R --no-restore --no-save
echo R script finished