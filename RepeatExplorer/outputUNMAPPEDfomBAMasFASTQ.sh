#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N exportUnmappedReadsFromBAM
#PBS -m abe
#--------------------------------------------------------------------------------------------
#Creates individual job files that output unmapped reads from BAM files in specified 'folder'
#and save the reads as pairs of 'fastq.gz' files
#These files can be later use to run RepeatExplorer
#Tomas Fer, 2026, tomas.fer@natur.cuni.cz
#v.0.0.2
#--------------------------------------------------------------------------------------------

#Specify folder (full path in MetaCentrum)
server=brno12-cerit #for output files
folder="/storage/${server}/home/${LOGNAME}/RE"

cd ${folder}
touch submitJobs.sh

#loop over all BAM files
for sample in $(ls *.bam | cut -d'.' -f1); do
	echo '#!/bin/bash' >> ${sample}.sh
	echo '#----------------MetaCentrum----------------' >> ${sample}.sh
	echo '#PBS -l walltime=04:00:00' >> ${sample}.sh
	echo '#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=16gb' >> ${sample}.sh
	echo '#PBS -j oe' >> ${sample}.sh
	echo '#PBS -o /storage/'"$server/home/$LOGNAME" >> ${sample}.sh
	echo '#PBS -N unmappedFromBAM_for_'"${sample}" >> ${sample}.sh
	echo 'folder='"${folder}" >> ${sample}.sh
	echo 'sample='"${sample}" >> ${sample}.sh
	echo 'cd ${SCRATCHDIR}' >> ${sample}.sh
	echo 'cp ${folder}/${sample}.bam .' >> ${sample}.sh
	echo 'module add samtools' >> ${sample}.sh
	echo 'module add fastq-pair' >> ${sample}.sh
	echo '#process BAM file' >> ${sample}.sh
	echo '#samtools view = -u: outputs uncompressed BAM, -f 12: select reads that have both Flag 4 (the read is unmapped) and Flag 8(the mate is unmapped) (4+8=12), -F 256: excludes secondary alignments' >> ${sample}.sh
	echo '#samtools collate = -u: expects uncompressed BAM input, -O: outputs BAM records that are properly paired and ordered, -: specifies reading from the standard input' >> ${sample}.sh
	echo '#samtools fastq = outputs correctly paired fastq.gz files' >> ${sample}.sh
	echo 'samtools view -u -f 12 -F 256 ${sample}.bam | samtools collate -u -O - | samtools fastq -1 ${sample}_unmapped_R1.fastq.gz -2 ${sample}_unmapped_R2.fastq.gz' >> ${sample}.sh
	echo '#decompress fastq.gz' >> ${sample}.sh
	echo 'pigz ${sample}_unmapped_R{1,2}.fastq.gz -d' >> ${sample}.sh
	echo '#check correct pairs using fastq_pair' >> ${sample}.sh
	echo 'fastq_pair ${sample}_unmapped_R1.fastq ${sample}_unmapped_R2.fastq' >> ${sample}.sh
	echo '#delete original fastq' >> ${sample}.sh
	echo 'rm ${sample}_unmapped_R{1,2}.fastq' >> ${sample}.sh
	echo '#delete files with single reads' >> ${sample}.sh
	echo 'rm ${sample}_unmapped_R{1,2}.fastq.single.fq' >> ${sample}.sh
	echo '#calculate nr. of reads' >> ${sample}.sh
	echo 'reads=$(expr $(cat ${sample}_unmapped_R1.fastq.paired.fq | wc -l) / 4)' >> ${sample}.sh
	echo 'echo -e "${sample}\t${reads}" > ${sample}_unmapped_nrReads.txt' >> ${sample}.sh
	echo '#compress correctly paired fastq files' >> ${sample}.sh
	echo 'pigz ${sample}_unmapped_R1.fastq.paired.fq -c > ${sample}_unmapped_R1.fastq.gz' >> ${sample}.sh
	echo 'pigz ${sample}_unmapped_R2.fastq.paired.fq -c > ${sample}_unmapped_R2.fastq.gz' >> ${sample}.sh
	echo '#delete uncompressed files' >> ${sample}.sh
	echo 'rm ${sample}_unmapped_R{1,2}.fastq.paired.fq' >> ${sample}.sh
	echo '#copy results home' >> ${sample}.sh
	echo 'mkdir ${folder}/${sample}' >> ${sample}.sh
	echo 'cp ${sample}_unmapped_R{1,2}.fastq.gz ${folder}/${sample}' >> ${sample}.sh
	echo 'cp ${sample}_unmapped_nrReads.txt ${folder}/${sample}' >> ${sample}.sh
	echo '#move BAM file' >> ${sample}.sh
	echo 'mv ${folder}/${sample}.bam ${folder}/${sample}' >> ${sample}.sh
	echo '#Clean scratch/work directory' >> ${sample}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> ${sample}.sh
	echo '  #delete scratch' >> ${sample}.sh
	echo '  if [[ ! $SCRATCHDIR == "" ]]; then' >> ${sample}.sh
	echo '    rm -rf $SCRATCHDIR/*' >> ${sample}.sh
	echo '  fi' >> ${sample}.sh
	echo 'fi' >> ${sample}.sh
	echo 'echo -e "\nFINISHED processing: ${sample}"' >> ${sample}.sh
	echo 'qsub '"${sample}"'.sh' >> submitJobs.sh
done

#make the submitter executable
chmod +x submitJobs.sh
echo -e "Go to '${folder}' and run 'submitJobs.sh'\n"
echo "Finished..."

