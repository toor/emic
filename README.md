# emic
Scripts to generate and analyse data used in the study: Entropic magnetic interlayer coupling (https://doi.org/10.1103/nfw6-jsn6)

Please note that data that you generate using this code will not be numerically identical to those presented in the paper, because the seed used for the random number generator may differ. However, the figures created using this data will be functionally identical.

## Prerequisites
To use this software, you need to install the Snakemake workflow management system (https://snakemake.readthedocs.io/en/stable/), which is most easily installed using the Anaconda package manager for Python (https://www.anaconda.com/docs/main). See the Anaconda and Snakemake documentation for details.
The data analysis scripts also require you to have the following Python packages installed:

* `numpy`
* `scipy`
* `matplotlib`
* `pandas`
* `numba`

Each of these packages can be installed into the `conda` environment at the same time as you install Snakemake.

The first step is to compile the binary:
`gcc -o main main.c`

Next, to run the Snakemake workflow, first edit `Snakefile` and replace `SIM_NAME` in the rules `all`, `run`, and `analyse_data`with the desired name of your output directory. Then, perform a Snakemake dry run using `snakemake -n`. This command tells you the number of jobs that Snakemake needs to run to generate these data.
To reproduce the data presented in the study, a large number of jobs need to be run, and we therefore recommend that you use a cluster software which can remotely run many jobs simultaneously. Snakemake integrates well with a variety of such software; for details, see: https://snakemake.readthedocs.io/en/stable/getting_started/migration.html and https://stackoverflow.com/questions/77929511/how-to-run-snakemake-8-on-a-slurm-cluster for use with Snakemake versions >= 8.0.

As it is setup at the moment, the Snakefile will generate data for small systems with a few interface spins, and plots the calculated entropy difference between parallel and antiparallel configurations of the macroscopic magnets for different values of $N_y$ as a function of $N_x$ (Figure 3 in the manuscript). `analyse_data.py` can also be modified to generate curves depicting magnetisation as a function of interfacial coupling strength $J_R$ for all combinations of $N_y$ and $N_x$.

For further information and/or clarification on the use of this software, please contact William Huddie (w.j.huddie@uu.nl)
