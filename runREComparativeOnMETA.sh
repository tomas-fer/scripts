#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=16:mem=128gb:scratch_local=64gb
#PBS -j oe
#PBS -N RepeatExplorerComparative
#PBS -m abe

#--------------------------------------------------------------------------------------------
#This script (executed on MetaCentrum via qsub) runs comparative RepeatExplorer analysis
#The list of samples is in the file 'listComparative.txt' in the $folder specified below,
#this file should contain three column
#1st column - the sample name
#2nd column - the sample abbreviation
#(all abbreviations have to be of same length and this length is specified as $prefixlength
#3rd column - genome size of the sample
#$suffix can be also specified
#Also specify read subsampling (number of random read pairs per sample)
#by setting the value of $subsample below
#
#Input files (pair of fastq.gz) - ${sample}${suffix}_R{1,2}.fastq.gz
#are in subfolders ($sample) in $folder
#i.e., this is easy to run in case the input was created with 'outputUNMAPPEDfomBAMasFASTQ.sh'
#
#After RepeatExplorer run, this script takes two output tables and uses 'plot_comparative_clustering_summary.R'
#to create a PDF summary plot (see https://github.com/kavonrtep/revis)
#(this requires 'optparse' package to be installed in /storage/$server/home/$LOGNAME/Rpackages44)
#Furthermore, the script extracts creates NJ trees for every cluster (based on observed/expected number of edges between species)
#see Vitales et al. 2020. Reconstructing phylogenetic relationships based on repeat sequence similarities
#and creates consensus network(s) using phangorn:::consensusNet() R function
#networks are created for 'prob' values from 0.05 to 0.15 (step 0.01)
#i.e., the proportion a split has to be present in all trees to be represented in the network
#this uses the R script 'REphylo.R' that has to be in /storage/home/${LOGNAME}/HybSeqSource
#(this requires 'phangorn' package to be installed in /storage/$server/home/$LOGNAME/Rpackages44)
#
#OUTPUT:
#- re_output_comparative.tar.gz        full results (compressed)
#- RE_comparative.log                  RepeatExplorer2 log
#- comparative_CLUSTER_TABLE.csv       automatic cluster annotation
#- COMPARATIVE_ANALYSIS_COUNTS.csv     nr. reads per cluster and sample
#- comparative_summary.png             basic summary plot
#- re_comparative_20colors.pdf         color summary plot reflecting sample genome size
#- OEnumberEdges                       folder with NJ trees and consensus network(s)
#
#Tomas Fer, 2025, tomas.fer@natur.cuni.cz
#v.0.0.2
#--------------------------------------------------------------------------------------------

#Specify this before running the script!!!
server=brno12-cerit
folder="/storage/brno12-cerit/home/tomasfer/Zingiberaceae_RE/BAM"
suffix="_unmapped"
subsample=100000
prefixlength=3

#Move to scratch
cd $SCRATCHDIR

export SINGULARITYENV_TMPDIR=$SCRATCHDIR
#download singularity image
echo "Downloading RE singularity image..."
singularity pull repex_tarean_0.3.12.sif library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e
#make singularity image
echo "Building RE singularity image..."
singularity build --sandbox repex_tarean repex_tarean_0.3.12.sif
module add seqtk/1.5

#prepare data (from filtered fastq.gz)
cp ${folder}/listComparative.txt .
cut -f1,2 listComparative.txt > list.txt
cut -f2,3 listComparative.txt > listGS.txt
#loop over samples
cat list.txt | while read -r sample prefix; do
	echo -en "Getting data..."
	echo ${sample}
	#copy data
	echo -e "\t...copy"
	cp ${folder}/${sample}/${sample}${suffix}_R{1,2}.fastq.gz .
	#sample to required coverage
	echo -e "\t...subsampling to ${subsample}"
	seqtk sample -s 10 ${sample}${suffix}_R1.fastq.gz ${subsample} > ${sample}${suffix}_R1.fastq
	seqtk sample -s 10 ${sample}${suffix}_R2.fastq.gz ${subsample} > ${sample}${suffix}_R2.fastq
	#merge PE fastq
	echo -e "\t...merge PE reads"
	seqtk mergepe ${sample}${suffix}_R1.fastq ${sample}${suffix}_R2.fastq > ${sample}${suffix}_merged.fastq
	#convert fastq to fasta
	echo -e "\t...convert to FASTA"
	seqtk seq -A ${sample}${suffix}_merged.fastq > ${sample}${suffix}_merged.fasta
	#add prefix
	echo -e "\t...add prefix...${prefix}\n"
	seqtk rename ${sample}${suffix}_merged.fasta ${prefix} > prefix_${sample}${suffix}_merged.fasta
	#clean
	rm ${sample}${suffix}_R{1,2}.fastq.gz ${sample}${suffix}_R{1,2}.fastq ${sample}${suffix}_merged.{fastq,fasta}
done
#concatenate
cat prefix* > comparative_final.fasta
rm prefix*
#run RepeatExplorer
echo "Running RE..."
singularity exec -e --bind ${PWD}:/data/ repex_tarean seqclust -C -l RE_comparative${subsample}.log -p --prefix_length ${prefixlength} -v /data/re_output_comparative${subsample} /data/comparative_final.fasta

#Make a plot using R (see https://github.com/kavonrtep/revis/tree/master)
echo -e "\nCreating plots..."
#cp ${folder}/plot_comparative_clustering_summary.R .
wget https://raw.githubusercontent.com/kavonrtep/revis/refs/heads/master/plot_comparative_clustering_summary.R
chmod +x plot_comparative_clustering_summary.R
#add R and path to library
export R_LIBS="/storage/$server/home/$LOGNAME/Rpackages44"
module add r/4.4.0-gcc-10.2.1-ssuwpvb
./plot_comparative_clustering_summary.R --cluster_table=re_output_comparative${subsample}/CLUSTER_TABLE.csv --comparative_counts=re_output_comparative${subsample}/COMPARATIVE_ANALYSIS_COUNTS.csv --number_of_colors=20 -g listGS.txt --output=re_comparative${subsample}_20colors.pdf

#extract data for phylogenetic tree and copy back home
nrsamples=$(cat list.txt | wc -l)
cd re_output_comparative${subsample}/seqclust/clustering/clusters
mkdir OEnumberEdges${subsample}
for i in $(ls -d dir*); do
	cp $i/observed_expected_number_of_edges.csv OEnumberEdges${subsample}/${i}.csv
done
cd OEnumberEdges${subsample}
#change acronymes to full names (in all *.csv files)
cat ../../../../../listComparative.txt | while read a b; do
	sed -i "s/$b/$a/" dir_CL*.csv
done
#make NJ trees for every cluster
cp /storage/${server}/home/${LOGNAME}/HybSeqSource/REphylo.R .
chmod +x REphylo.R
./REphylo.R $nrsamples
cd ..
cp -r OEnumberEdges${subsample} ${folder}/
cd ../../../..

#copy comparative and annotation results back home
cp RE_comparative${subsample}.log ${folder}/
cp re_output_comparative${subsample}/CLUSTER_TABLE.csv ${folder}/comparative_${subsample}_CLUSTER_TABLE.csv
cp re_output_comparative${subsample}/COMPARATIVE_ANALYSIS_COUNTS.csv ${folder}/COMPARATIVE_${subsample}_ANALYSIS_COUNTS.csv
cp re_output_comparative${subsample}/comparative_summary.png ${folder}/comparative_summary${subsample}.png
cp re_comparative${subsample}_20colors.pdf ${folder}/

#gzip the result directory
echo "Gzip result folder..."
tar -czf re_output_comparative${subsample}.tar.gz re_output_comparative${subsample}

#copy back home
echo "Copying data home..."
cp re_output_comparative${subsample}.tar.gz ${folder}/

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
fi
exit
