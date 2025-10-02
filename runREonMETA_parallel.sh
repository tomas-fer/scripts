#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N RepeatExplorer
#PBS -m abe
#--------------------------------------------------------------------------------------------
#Creates individual job files for running RepeatExplorer
#The list of samples is in the file 'list.txt' in the $folder specified below,
#suffix can be also specified
#Input files (pair of fastq.gz) - ${sample}${suffix}_R{1,2}.fastq.gz
#are in subfolders ($sample) in $folder
#i.e., this is easy to run in case the input was created with 'outputUNMAPPEDfomBAMasFASTQ.sh'
#OUTPUT:
#- ${sample}.tar.gz
#- RE_${sample}.log
#- ${sample}_CLUSTER_TABLE.csv
#
#Tomas Fer, 2025, tomas.fer@natur.cuni.cz
#v.0.0.1
#--------------------------------------------------------------------------------------------

server=brno2
folder="/storage/brno12-cerit/home/tomasfer/Zingiberaceae_RE/BAM"
suffix="_unmapped"
cd ${folder}
touch submitREjobs.sh
#loop over list of samples
for sample in $(cat list.txt); do
	echo '#!/bin/bash' >> RE_${sample}.sh
	echo '#----------------MetaCentrum----------------' >> RE_${sample}.sh
	echo '#PBS -l walltime=24:00:00' >> RE_${sample}.sh
	echo '#PBS -l select=1:ncpus=16:mem=128gb:scratch_local=64gb' >> RE_${sample}.sh
	echo '#PBS -j oe' >> RE_${sample}.sh
	echo '#PBS -o /storage/'"$server/home/$LOGNAME" >> RE_${sample}.sh
	echo '#PBS -N RE_for_'"${sample}" >> RE_${sample}.sh
	echo 'export SINGULARITYENV_TMPDIR=$SCRATCHDIR' >> RE_${sample}.sh
	echo 'cd $SCRATCHDIR' >> RE_${sample}.sh
	echo 'folder='"${folder}" >> RE_${sample}.sh
	echo 'sample='"${sample}" >> RE_${sample}.sh
	echo 'suffix='"${suffix}" >> RE_${sample}.sh
	echo '#download singularity image' >> RE_${sample}.sh
	echo 'echo -e "\nDownloading RE singularity image..."' >> RE_${sample}.sh
	echo 'singularity pull repex_tarean_0.3.12.sif library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e' >> RE_${sample}.sh
	echo '#make singularity image' >> RE_${sample}.sh
	echo 'echo -e "\nBuilding RE singularity image..."' >> RE_${sample}.sh
	echo 'singularity build --sandbox repex_tarean repex_tarean_0.3.12.sif' >> RE_${sample}.sh
	echo '#add module(s)' >> RE_${sample}.sh
	echo 'module add seqtk/1.5' >> RE_${sample}.sh
	echo '#prepare data (from filtered fastq.gz)' >> RE_${sample}.sh
	echo 'echo' >> RE_${sample}.sh
	echo 'echo -en "Getting data..."' >> RE_${sample}.sh
	echo 'echo ${sample}' >> RE_${sample}.sh
	echo '#copy data' >> RE_${sample}.sh
	echo 'echo -e "\t...copy"' >> RE_${sample}.sh
	echo 'cp ${folder}/${sample}/${sample}${suffix}_R{1,2}.fastq.gz .' >> RE_${sample}.sh
	echo '#unzip fastq.gz' >> RE_${sample}.sh
	echo 'echo -e "\t...unzip"' >> RE_${sample}.sh
	echo 'gzip -d ${sample}${suffix}_R{1,2}.fastq.gz' >> RE_${sample}.sh
	echo '#merge PE fastq' >> RE_${sample}.sh
	echo 'echo -e "\t...merge PE reads"' >> RE_${sample}.sh
	echo 'seqtk mergepe ${sample}${suffix}_R1.fastq ${sample}${suffix}_R2.fastq > ${sample}${suffix}_merged.fastq' >> RE_${sample}.sh
	echo '#convert fastq to fasta' >> RE_${sample}.sh
	echo 'echo -e "\t...convert to FASTA"' >> RE_${sample}.sh
	echo 'seqtk seq -A ${sample}${suffix}_merged.fastq > ${sample}${suffix}_merged.fasta' >> RE_${sample}.sh
	echo '#run RepeatExplorer (and create a logfile)' >> RE_${sample}.sh
	echo 'echo -e "\nRunning RE..."' >> RE_${sample}.sh
	echo 'singularity exec -e --bind ${PWD}:/data/ repex_tarean seqclust -C -l RE_${sample}.log -p -v /data/${sample} /data/${sample}${suffix}_merged.fasta' >> RE_${sample}.sh
	echo '#copy annotation history and logfile back home' >> RE_${sample}.sh
	echo 'cp ${sample}/CLUSTER_TABLE.csv ${folder}/${sample}/${sample}_CLUSTER_TABLE.csv' >> RE_${sample}.sh
	echo 'cp RE_${sample}.log ${folder}/${sample}' >> RE_${sample}.sh
	echo '#gzip the result directory' >> RE_${sample}.sh
	echo 'echo -e "\nGzip result folder..."' >> RE_${sample}.sh
	echo 'tar -czf ${sample}.tar.gz ${sample}' >> RE_${sample}.sh
	echo '#copy back home' >> RE_${sample}.sh
	echo 'echo -e "\nCopying data home..."' >> RE_${sample}.sh
	echo 'cp ${sample}.tar.gz ${folder}/${sample}' >> RE_${sample}.sh
	echo '#Clean scratch/work directory' >> RE_${sample}.sh
	echo 'if [[ $PBS_O_HOST == *".cz" ]]; then' >> RE_${sample}.sh
	echo '  #delete scratch' >> RE_${sample}.sh
	echo '  if [[ ! $SCRATCHDIR == "" ]]; then' >> RE_${sample}.sh
	echo '    rm -rf $SCRATCHDIR/*' >> RE_${sample}.sh
	echo '  fi' >> RE_${sample}.sh
	echo 'fi' >> RE_${sample}.sh
	echo 'echo -e "\nFINISHED processing: ${sample}"' >> RE_${sample}.sh
	echo 'qsub RE_'"${sample}"'.sh' >> submitREjobs.sh
done

#make the submitter executable
chmod +x submitREjobs.sh
echo -e "Go to '${folder}' and run 'submitREjobs.sh'\n"
echo finished
