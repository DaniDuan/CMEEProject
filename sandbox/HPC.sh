#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#module load anaconda3/personal
echo "running HPC"
python ./run_on_HPC.py
echo "finished running"
# this is a comment at the end of the file
