#!/bin/bash
#$ -N TypeI100           # <- name your job
#$ -j y              # <- join output and error for easy reading
#$ -pe openmpi 30     # <- requesting 16 cores
#$ -l m_mem_free=2G  # <- start with MUCH LESS RAM

time Rscript --no-save BodyTypeI.R
