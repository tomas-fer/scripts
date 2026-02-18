#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N RepeatExplorer
#PBS -m abe
#-------------------------------------------
#Summarizes number of unmapped reads from BAM files

folder="/storage/brno12-cerit/home/tomasfer/PicrisOlympica_RE"

cd ${folder}
for i in $(ls -l | grep ^d | awk '{ print $9 }'); do
	cat ${i}/*.txt >> unmappedReadsSummary.txt
done
