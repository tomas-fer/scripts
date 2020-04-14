#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N cpDNAmapping
#PBS -m abe

#-------------------HYDRA-------------------
#$ -S /bin/bash
#$ -q sThC.q
#$ -l mres=1G
#$ -cwd
#$ -j y
#$ -N cpDNAmapping
#$ -o cpDNAmapping.log

#------------------------------------------------------------------------------------
#Tomas Fer, 2018, 2019
#Department of Botany, Charles University, Prague, Czech Republic
#tomas.fer@natur.cuni.cz
#------------------------------------------------------------------------------------

#Maps reads to the full plastome reference (integrated within HybPhyloMaker)
#(1) BWA mapping of filtered reads to cpDNA reference (using filtered reads without duplicates)
#(2) Consensus call with kindel
#(3) Combine sequences to single FASTA file
#(4) Calculate proportion of missing data per samples

#USAGE: ./cpDNA_mappingMETA.sh (or qsub cpDNA_mappingMETA.sh)

#Requirements:
#BWA
#samtools
#kindel v.0.1.4 (pip3 install 'kindel==0.1.4')

if [[ $PBS_O_HOST == *".cz" ]]; then
	echo -e "\ncpDNAmapping is running on MetaCentrum..."
	#settings for MetaCentrum
	#Move to scratch
	cd $SCRATCHDIR
	#Copy file with settings from home and set variables from settings.cfg
	cp $PBS_O_WORKDIR/settings.cfg .
	. settings.cfg
	. /packages/run/modules-2.0/init/bash
	path=/storage/$server/home/$LOGNAME/$data
	source=/storage/$server/home/$LOGNAME/HybSeqSource
	#Add necessary modules
	module add bwa-0.7.15
	module add samtools-1.3
	module add python34-modules-gcc #adds also kindel
elif [[ $HOSTNAME == compute-*-*.local ]]; then
	echo -e "\ncpDNAmapping is running on Hydra..."
	#settings for Hydra
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir cpDNAmapping
	cd cpDNAmapping
	#Add necessary modules
	module load bioinformatics/bwa/0.7.12
	module load bioinformatics/samtools/1.3
	module load bioinformatics/anaconda3/5.1 #adds also kindel
else
	echo -e "\nHybPhyloMaker9 is running locally..."
	#settings for local run
	#set variables from settings.cfg
	. settings.cfg
	path=../$data
	source=../HybSeqSource
	#Make and enter work directory
	mkdir cpDNAmapping
	cd cpDNAmapping
fi

#Copy data (all deduplicated FASTQ files from from all samples)
find $path/20filtered/ -name '*no-dups*' -type f -exec cp '{}' . \;
#Copy reference
cp $source/$cpDNA .

#Make folder for results
mkdir $path/fullplastome

#Index reference
bwa index $cpDNA #index reference
ls *-1P_no-dups* | cut -d"-" -f1,2 > list.txt
for file in $(cat list.txt); do
	echo -e "Processing $file"
	echo -e "...mapping"
	bwa mem $cpDNA ${file}-1P_no-dups.fastq.gz ${file}-2P_no-dups.fastq.gz > ${file}.sam #mapping
	echo -e "...converting"
	samtools view -bS -o ${file}.bam ${file}.sam #convert to bam
	nrmapped=$(samtools view -F 0x904 -c ${file}.bam) #number of mapped reads in BAM
	nrall=$(samtools view -c ${file}.bam) #number of all reads in BAM
	percmapped=$(echo -e "scale=5;100 * ($nrmapped / $nrall)" | bc) #calculate percentage of mapped reads
	echo -e "$file\t$nrall\t$nrmapped\t$percmapped" >> mapping_summary.txt
	cp mapping_summary.txt $path/fullplastome
	echo -e "...sorting"
	samtools sort ${file}.bam -o ${file}_sorted.bam
	echo -e "...indexing"
	samtools index ${file}_sorted.bam
	cp ${file}_sorted.bam $path/fullplastome/${file}.bam
	cp ${file}_sorted.bam.bai $path/fullplastome/${file}.bam.bai
	kindel -m $mincov -t $majthres ${file}_sorted.bam > ${file}.fasta
	#change name in fasta file
	sed -i '1d' ${file}.fasta #delete first line
	echo ">$file" > header.txt
	cat header.txt ${file}.fasta > tmp && mv tmp ${file}.fasta
	rm header.txt ${file}.sam ${file}.bam #cleaning
	#Remove line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${file}.fasta > tmp && mv tmp ${file}.fasta
	cp ${file}.fasta $path/fullplastome/${file}.fasta
done
#combine fasta files
cat *.fasta > consensus.fasta
#Remove line breaks from fasta file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' consensus.fasta > tmp && mv tmp consensus.fasta
cp consensus.fasta $path/fullplastome/fullplastomes_consensus.fasta

#Calculate percentage of missing data per accession
#Calculate length of alignment: 1. get second line and count length, 2. decrease value by one (because previous command also counted LF)
length=$(cat consensus.fasta | head -n 2 | tail -n 1 | wc -c)
length=`expr $length - 1`
#Replace newline with ' ' if line starts with '>' (i.e., merge headers with data into single line separated by space)
cat consensus.fasta | sed '/^>/{N; s/\n/ /;}' > consensus.modif.fasta
#Cut first part until space, i.e. header, and remove '>'
cat consensus.modif.fasta | cut -f1 -d" " | sed 's/>//' > headers.txt
#Cut only part after the first space, i.e., only sequence, change all missing data (-, ?, N) to 'n', replace all other characters then 'n' by nothing and print percentage of 'n's in each sequence
cat consensus.modif.fasta | cut -f2 -d" " | sed 's/[?N-]/n/g' | sed 's/[^n]//g' | awk -v val=$length '{ print (length*100)/val }' > missingpercentage.txt
paste headers.txt missingpercentage.txt > fullplastomes_missingperc.txt
#Calculate mean of all values
echo -e "MEAN\t$(awk '{ sum += $2; n++ } END { if (n > 0) print sum / n; }' fullplastomes_missingperc.txt)" > mean.txt
cat fullplastomes_missingperc.txt mean.txt > tmp && mv tmp fullplastomes_missingperc.txt
rm consensus.modif.fasta headers.txt missingpercentage.txt mean.txt
#Copy results to home
cp fullplastomes_missingperc.txt $path/fullplastome

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
else
	cd ..
	rm -r cpDNAmapping
fi

echo -e "\ncpDNA mapping finished...\n"
