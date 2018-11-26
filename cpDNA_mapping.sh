#(1) BWA mapping of filtered reads to cpDNA reference (use filtered reads without duplicates from HybPhyloMaker, copy them to a single folder)
#(2) Consensus call with kindel
#(3) Combine sequences to single FASTA file

#Requirements:
#BWA
#samtools
#kindel v.0.1.4 (pip3 install 'kindel==0.1.4')

#USAGE: ./cpDNA_mapping.sh

#------------------------------------------------------------------------------------
#Tomas Fer, 2018
#Department of Botany, Charles University, Prague, Czech Republic
#tomas.fer@natur.cuni.cz
#------------------------------------------------------------------------------------

#SETTINGS
reference=cpDNA.fas
mincov=2
majthres=0.9

#Index reference
bwa index $reference #index reference
ls *-1P_no-dups* | cut -d"-" -f1,2 > list.txt
for file in $(cat list.txt); do
	echo -e "Processing $file"
	echo -e "...mapping"
	bwa mem $reference ${file}-1P_no-dups.fastq.gz ${file}-2P_no-dups.fastq.gz > ${file}.sam #mapping
	echo -e "...converting"
	samtools view -bS -o ${file}.bam ${file}.sam #convert to bam
	nrmapped=$(samtools view -F 0x04 -c ${file}.bam) #number of mapped reads in BAM
	nrall=$(samtools view -c ${file}.bam) #number of all reads in BAM
	percmapped=$(echo -e "scale=5;100 * ($nrmapped / $nrall)" | bc) #calculate percentage of mapped reads
	echo -e "$file\t$nrall\t$nrmapped\t$percmapped" >> mapping_summary.txt
	echo -e "...sorting"
	samtools sort ${file}.bam ${file}_sorted
	echo -e "...indexing"
	samtools index ${file}_sorted.bam
	kindel -m $mincov -t $majthres ${file}_sorted.bam > ${file}.fasta
	#change name in fasta file
	sed -i '1d' ${file}.fasta #delete first line
	echo ">$file" > header.txt
	cat header.txt ${file}.fasta > tmp && mv tmp ${file}.fasta
	rm header.txt ${file}.sam ${file}.bam #cleaning
	#Remove line breaks from fasta file
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${file}.fasta > tmp && mv tmp ${file}.fasta
done
#combine fasta files
cat *.fasta > consensus.fasta
#Remove line breaks from fasta file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' consensus.fasta > tmp && mv tmp consensus.fasta
