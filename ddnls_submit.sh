#!/bin/bash
#SBATCH -D /scratch/users/ccastrocastr/ddnls_abc4
#SBATCH -p serial
#SBATCH -N 10
#SBATCH --mail-user ccastrocastr@smu.edu 
#SBATCH --mail-type=all

# run the code
time ./driver.exe > diary_abc4.txt
gprof driver.exe > profiling_data.txt

module load MATLAB
matlab plots_abc4
