#!/bin/bash

#SBATCH --ntasks=16
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --job-name=PE_chip_snakemake
#SBATCH --account=PES0738

module load python/3.7-2019.10

source activate snakemake

cd /fs/project/PES0738/testing/chip_snakemake

snakemake --cores 16 -s chipseq_bowtie2_pe.Snakefile --use-conda --conda-frontend conda --nolock --rerun-incomplete 