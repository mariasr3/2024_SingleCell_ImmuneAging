#!/bin/bash
#SBATCH --chdir=/gpfs/projects/bsc83/Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/14_SenCID
#SBATCH --cpus-per-task=64
#SBATCH --qos=gp_debug
#SBATCH --time=02:00:00
#SBATCH --job-name=14.SenCID
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --account=bsc83
#SBATCH --constraint=highmem
#SBATCH --output=/gpfs/projects/bsc83/Projects/scRNAseq/msopena/02_OneK1K_Age/scripts/Analysis/log_dir/outputs/14.SenCID.%j.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/scRNAseq/msopena/02_OneK1K_Age/scripts/Analysis/log_dir/errors/14.SenCID.%j.err

module purge
module load intel impi/ bsc/ mkl/2024.0 geos/3.12.1 hdf5/ proj/ gcc/ gsl R/4.3.2  transfer/1.0 szip miniconda
conda init
conda activate SenCID
SenCID --filepath /gpfs/projects/bsc83/Projects/scRNAseq/msopena/02_OneK1K_Age/robjects/14_SenCID/ --denoising t --fileclass txt --binarize t --save_tmps save_tmps
#SenCID --filepath . --denoising t --fileclass txt --binarize t --save_tmps save_tmps