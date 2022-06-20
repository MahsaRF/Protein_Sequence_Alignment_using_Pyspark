#! /bin/bash
# Mahsa Rezaei firuzkuhi
#SBATCH -J prot
#SBATCH -t 01:00:00
#SBATCH -N 5
#SBATCH -p cosc6339
#SBATCH --cpus-per-task 2

 
cd /home2/stud25/
module load spark/2.3.4
crill-spark-submit -a 2 proteinscore-with-text.py protein-1.txt db-large.txt 2 output
exit
