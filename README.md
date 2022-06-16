[![dry-run](https://github.com/admixVIE/sstar-analysis/actions/workflows/dry-run.yml/badge.svg?branch=main)](https://github.com/admixVIE/sstar-analysis/actions)

# sstar-analysis

## Introduction

This repo contains [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines for replicating the analysis in our `sstar` paper. These pipelines were tested on Linux operating systems (CentOS8, Oracle Linux 8, and Ubuntu 20.04).

To replicate our analysis, users should install [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) first, then use the following commands to create a virtual environment for the analysis.

	conda config --set safety_checks disabled
	conda config --set channel_priority strict
	conda env create -f environment.yml
	conda activate sstar-analysis

For tools that cannot be installed through `conda`, users could follow the commands below.

	mkdir ext && cd ext

	# Download SPrime and its pipeline
	mkdir SPrime && cd SPrime
	wget https://faculty.washington.edu/browning/sprime.jar
	git clone https://github.com/YingZhou001/sprimepipeline
	chmod a+x sprimepipeline/pub.pipeline.pbs/tools/map_arch_genome/map_arch
	sed 's/out<-c()/out<-data.frame()/' sprimepipeline/pub.pipeline.pbs/tools/score_summary.r > tmp
	mv tmp sprimepipeline/pub.pipeline.pbs/tools/score_summary.r
	cd ..

	# Download SkovHMM
	git clone https://github.com/LauritsSkov/Introgression-detection SkovHMM

	# Download ArchaicSeeker2.0
	git clone https://github.com/Shuhua-Group/ArchaicSeeker2.0
	cd ArchaicSeeker2.0
	cd libnlopt.so.0_free
	chmod a+x ArchaicSeeker2
	cd ../..

Users also need to download `ms.tar.gz` from [Hudson Lab](http://home.uchicago.edu/~rhudson1/source/mksamples.html) and decompress it under the `ext` folder and compile it with the following commands.

	cd msdir
	${CONDA_PREFIX}/bin/gcc -o ms ms.c streec.c rand1.c -lm
	cd ../..

## Running the pipelines

After installing the tools above, users could test the pipelines locally by using dry-run.

	snakemake -np

If users want to run all the pipelines, users can use the following command.

	snakemake -c 1 --use-conda

`-c` specifies the number of threads and `snakemake` could run jobs parallelly as many as possible with the given number of threads.

However, we recommend users to run each pipeline individually. Users need to run the two pipelines for simulation first.

	snakemake -s workflows/1src/simulation.snake -c 1
	snakemake -s workflows/2src/simulation.snake -c 1

Other pipelines, including `workflows/1src/sstar.snake`, `workflows/1src/sprime.snake`, `workflows/1src/skovhmm.snake`, `workflows/2src/sstar.snake`, `workflows/2src/sprime.snake`, `workflows/2src/archaicseeker2.snake`, can be executed in any order after simulation.

For the SkovHMM pipeline, users need to add the `--use-conda` argument, because it requires `Python2`, which is different from the main environment `sstar-analysis`. A specific environment for SkovHMM needs to be created.

	snakemake -s workflows/1src/skovhmm.snake -c 1 --use-conda

Finally, users could plot the results.

	snakemake -s workflows/plot.snake -c 1

The plots and tables are in `results/plots`. They may be slightly different from those in our manuscript because of random effects.

Users could also submit the pipelines to HPC. Users should create [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) depending on the HPC scheduler. Users could find an example profile for `SLURM` in `config/slurm/config.yaml`. To submit jobs by `SLURM`, users first create a folder `logs_slurm`.

	mkdir logs_slurm
	
Then, for example, to run simulation in the cluster,

	snakemake -s workflows/1src/simulation.snake --profile config/slurm -j 200 

`-j` specifies the number of threads in the cluster.
