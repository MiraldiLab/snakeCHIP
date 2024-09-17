#!/bin/bash

#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --job-name=SE_chip_snakemake
#SBATCH --account=PES0738

module load python/3.7-2019.10

source activate snakeatac

cd /fs/ess/PES0738/20220614_maxatac_v1_data/snakemake/chip

snakemake --cores 14 -s chipseq_bowtie2_se.Snakefile --use-conda --conda-frontend conda --configfile ./inputs/config.yaml --nolock --rerun-incomplete --latency-wait 30