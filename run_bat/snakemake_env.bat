#!/bin/bash

#SBATCH --ntasks=16
#SBATCH --time=48:00:00
#SBATCH --job-name=SE_chip_snakemake
#SBATCH --account=PES0738

module load miniconda3/24.1.2-py310

source activate snakemake

cd /fs/ess/PES0738/20220614_maxatac_v1_data/snakemake/chip

snakemake --cores 16 --use-conda --conda-frontend conda --conda-create-envs-only --configfile ./inputs/config.yaml -s chipseq_bowtie2_se.Snakefile 
