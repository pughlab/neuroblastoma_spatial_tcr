#!/bin/bash
#SBATCH -J jobname
#SBATCH -t 5-00:00:00
#SBATCH -p himem
#SBATCH --mem 50G
#SBATCH -o %x.out
#SBATCH -c 1
#SBATCH -N 1

module load java/18
module load mixcr/4.6.0

mixcr -Xmx50g analyze local:cap-tcr /path/to/rawdata_R1.fastq.gz /path/to/rawdata_R2.fastq.gz /path/to/mixcr4/output/outfile
