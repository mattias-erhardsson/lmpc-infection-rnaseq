# lmpc-infection-rnaseq
Supplemental files, other data and analysis for the manuscript "The mouse gastric surface epithelial cell and its response to early Helicobacter pylori infection".

The supplemental files are in the supplemental files folder https://github.com/mattias-erhardsson/lmpc-infection-rnaseq/tree/main/Virulence%20supplemental%20files.

R version 4.3.2 is used. Renv is used for package version control. In script 4 you need to launch Cytoscape before running through the script. Script 7 is a jupyter notebook python script, not an R script.

Some results in the output folders might vary from the results presented in the article. If this is the case then the reason is that the script has been run at different times. The files in the commits 2025-05-08 are from a script re-run on the the first submission date on a system with freeshly installed R 4.3.2 with no packages installed, and with commit c2007e2 freshly cloned, in order to ensure other people can run the script and to help make the analysis as reproducible as possible. The reason why results might vary between runs is that there may be dependency packages and linked databases not fully defined in this github repository that may have been updated between script runs. In order to more accurately replicate results from the article, consider replacing input generated throughout the script with what is provided in the supplemental files for the manuscript.
