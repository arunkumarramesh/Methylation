#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100000
#SBATCH -J rename_pe
#SBATCH -o rename_pe.out
#SBATCH -e rename_pe.err

cd /proj/popgen/a.ramesh/projects/methylomes/arabidopsis/data/

cat SRR3300146_1.paired.fq.gz SRR3300147_1.paired.fq.gz > 7067_1.paired.fq.gz
cat SRR3300312_1.paired.fq.gz SRR3300313_1.paired.fq.gz SRR3300314_1.paired.fq.gz > 7521_1.paired.fq.gz
cat SRR3300342_1.paired.fq.gz SRR3300343_1.paired.fq.gz > 8239_1.paired.fq.gz
cat SRR3300146_2.paired.fq.gz SRR3300147_2.paired.fq.gz > 7067_2.paired.fq.gz
cat SRR3300312_2.paired.fq.gz SRR3300313_2.paired.fq.gz SRR3300314_2.paired.fq.gz > 7521_2.paired.fq.gz
cat SRR3300342_2.paired.fq.gz SRR3300343_2.paired.fq.gz > 8239_2.paired.fq.gz
