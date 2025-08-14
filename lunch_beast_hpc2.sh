#!/bin/bash
#SBATCH --job-name=lunch_beast_hpc2
#SBATCH --ntasks=1
#SBATCH -p gdec
#SBATCH --time=10-00:00:00
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=XXXXX
#SBATCH --mail-type=ALL 
#SBATCH --error="XXXXX/lunch_beast_hpc2.err"
#SBATCH --output="XXXXX/lunch_beast_hpc2.out"

# To download the correct version of Beast (2.7.8) (it is important not to take another version to avoid conflicts):
# https://www.beast2.org

# If linux:
# wget https://www.beast2.org/download-linux-x86/
# tar fxz BEAST.v2.7.7.Linux.x86.tgz
# cd beast/bin
# then run beast: ./beast/bin/beast

module load conda/4.12.0
source ~/.bashrc
module load java/oracle-1.8.0_45
module load gcc/8.1.0
module load beagle/3.1.2

WORKING_DIRECTORY=XXXXX/beast/bin

cd $WORKING_DIRECTORY

./beast XXXXX/Murat_parameters_cleaned_calibrated_all_prior.xml
