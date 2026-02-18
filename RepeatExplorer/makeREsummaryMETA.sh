#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=1:0:0
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -j oe
#PBS -N RE_summary
#PBS -m abe
#--------------------------------------------------------------------------------------------
#Summarize results of RepeatExplorer2 runs for multiple samples in specified '$folder'
#Calculate number of reads/proportion of all identified repeat types in RepeatExplorer cluster annotation file
#Based on annotations in the 5th (automatic) or 7th (final) column
#Annotation must follow the standards of RepeatExplorer, i.e. always start with 'All' etc.
#Proportion is calculated based on number of nuclear reads (i.e. after cp/mt reads removal)
#Size is corrected to the genome size in Mb/1C (provide an integer, i.e. whole number)
#Requires a file 'listGS.txt' with two TAB separated columns: sample name and genome size (1C in Mbp)
#(stored in '$folder'
#
#Tomas Fer, 2022-2026
#tomas.fer@natur.cuni.cz
#Version: 0.6 (2026-02-18)
#--------------------------------------------------------------------------------------------

#Local usage: ./makeREsummary.sh annotationtype filename genomeSizeMb1C
#e.g. ./makeREsummary.sh a CLUSTER_TABLE.csv 1361 #summary of automatic annotation (i.e., 5th column)
#e.g. ./makeREsummary.sh f CLUSTER_TABLE.csv 1361 #summary of final annotation (i.e., 7th column)

#Output files:
#filename_{automatic,final}_annotation.txt - nr. reads and annotation of all top clusters (columns 4 and 5/7 from filename)
#filename_stat.txt - read statistics
#filename_sumprop.txt - summary of reads/proportions/sizes(Mb) of all repeat types
#filename_sumpropsimple.txt - the same as the previous but with simplified repeat names
#filename_sumpropgen.txt - summary of reads/proportions/sizes(Mb) of basic repeat types
#RE_stat.txt - read statistics for all samples
#RE_sumpropgenMb.txt - summary of basic repeat types for all samples
#RE_sumpropsimpleMb.txt - summary of all repeat types for all samples
#REgen.pdf - barplot of basic repeat types for all samples
#REgenPerc.pdf - percentage barplot of basic repeat types for all samples
#REsimple.pdf - barplot of all repeat types for all samples
#REsimplePerc.pdf - percentage barplot of all repeat types for all samples

#Specify folder (full path in MetaCentrum)
server=brno12-cerit
folder="/storage/${server}/home/${LOGNAME}/RE"
annot=a

cd ${folder}

#Check R script
if [ ! -f "REplot.R" ]; then
	echo -e "Downloading REplot.R\n"
	wget https://raw.githubusercontent.com/tomas-fer/scripts/refs/heads/master/RepeatExplorer/REplot.R
	echo 
fi

#loop over a list of samples
cat listGS.txt | while read -r sample gs; do
	cd $sample
	#Check if input file exists
	if [ ! -f "${sample}_CLUSTER_TABLE.csv" ]; then
		echo -e "Error: input file for $sample does not exist" >&2; cd ..; continue
	fi
	#Check if genome size is an integer
	re='^[0-9]+$'
	if ! [[ $gs =~ $re ]] ; then
		echo -e "error: GS of $sample is not an integer" >&2; cd ..; continue
	fi
	#Check if annotation type is 'f' or 'a'
	if ! [[ $annot =~ ^f$ || $annot =~ ^a$ ]] ; then
		echo "error: annotation type should be 'f' or 'a'" >&2; cd ..; exit 3
	fi
	#Get appropriate annotation column
	if [[ $annot =~ ^a$ ]] ; then
		annotation=automatic
		#get 4th and 5th columns only, i.e., read numbers and automatic annotation; remove lines with TABs only
		awk '{ print $4"\t"$5"\t" }' ${sample}_CLUSTER_TABLE.csv | sed '/^[[:blank:]]/d' > annotation.txt
	elif [[ $annot =~ ^f$ ]] ; then
		annotation=final
		#get 4th and 7th columns only, i.e., read numbers and final annotation; remove lines with TABs only
		awk -F'\t' '{ print $4"\t"$7"\t" }' ${sample}_CLUSTER_TABLE.csv | sed '/^[[:blank:]]/d' > annotation.txt
	fi
	#Check if annotations exist
	if [ -z "$(awk '{ print $2 }' annotation.txt | grep "All")" ]; then
		echo "error: annotation for ${sample} missing or not in the correct format" >&2; cd ..; continue
	fi
	#Initial screen output
	echo
	echo "File: ${sample}_CLUSTER_TABLE.csv"
	echo "Annotation: $annotation"
	echo "Genome size (1C): $gs Mb"
	echo
	#Remove path (i.e. everything before the last '/') and then the suffix from file name
	#file=$(sed 's@.*/@@' <<< ${sample}_CLUSTER_TABLE.csv | cut -d'.' -f1)
	file=$sample
	#Get unique repeat annotations
	echo "Summarizing repeat types..."
	awk '{print $2 }' annotation.txt | sort | uniq | grep "All" > types.txt
	for type in $(cat types.txt); do
		#echo ${type}
		echo -en "${type}\t" >> ${file}_summary.txt
		#Summarize number of reads on lines matching specific repeat
		grep -P "[0-9]\t${type}\t" annotation.txt | awk 'BEGIN{FS="\t"; sum=0} {sum+=$1} END{print sum}' >> ${file}_summary.txt
	done
	#Read statistics
	echo "Repeat statistics..."
	all=$(grep "Number_of_analyzed_reads" ${sample}_CLUSTER_TABLE.csv | awk '{ print $2 }')
	plastid=$(grep "All/organelle/plastid" ${file}_summary.txt | awk '{ print $2 }')
	[ -z "$plastid" ] && plastid=0 #set 'plastid' to '0' if empty
	mitochondrial=$(grep "All/organelle/mitochondria" ${file}_summary.txt | awk '{ print $2 }')
	[ -z "$mitochondrial" ] && mitochondrial=0 #set 'mitochondrial' to '0' if empty
	allnuclear=$(( $all - $plastid - $mitochondrial ))
	singlets=$(grep "Number_of_singlets" ${sample}_CLUSTER_TABLE.csv | awk '{ print $2 }')
	topclust=$(grep -P "All" annotation.txt | awk 'BEGIN{FS="\t"; sum=0} {sum+=$1} END{print sum}')
	topclustnuc=$(( $topclust - $plastid - $mitochondrial ))
	minorclust=$(( $all - $topclust - $singlets))
	#Remove plastid and mitochondria from types.txt
	grep -vE "plastid|mitochondria" types.txt > tmp && mv tmp types.txt
	grep -vE "plastid|mitochondria" ${file}_summary.txt > ${file}_summary2.txt
	#Add small clusters and singlets
	echo -e "smallClusters\t${minorclust}" >> ${file}_summary2.txt
	echo -e "singlets\t${singlets}" >> ${file}_summary2.txt
	echo -e "smallClusters\nsinglets" >> types.txt
	#Calculate read proportion and length in Mb
	echo "Calculating read proportion and size in Mb..."
	for type in $(cat types.txt); do
		reads=$(grep -P "${type}\t" ${file}_summary2.txt | awk '{ print $2 }')
		perc=$(echo -e "scale=8;100 * ($reads / $allnuclear)" | bc)
		mb=$(echo -e "scale=8;$gs * ($perc / 100)" | bc)
		echo -e "${perc}\t${mb}" >> ${file}_summary3.txt
	done
	mv annotation.txt ${file}_${annotation}_annotation.txt
	#Combine tables
	paste ${file}_summary2.txt ${file}_summary3.txt >> ${file}_sumprop.txt
	rm ${file}_summary.txt ${file}_summary2.txt ${file}_summary3.txt types.txt
	#Table with simplified repeat names
	cat ${file}_sumprop.txt | cut -f1 | sed 's@.*/@@' > 1.txt #removes everything before the last '/'
	cat ${file}_sumprop.txt | cut -f2-4 > 2.txt
	paste 1.txt 2.txt > ${file}_sumpropsimple.txt
	rm 1.txt 2.txt
	#Remove double quotations
	sed -i 's/"//g' ${file}_sumprop.txt
	sed -i 's/"//g' ${file}_sumpropsimple.txt
	#General summary table
	echo "Summary tables..."
	sed 's/All\/repeat/All/' ${file}_sumprop.txt > ${file}_sumprop2.txt #change 'All/repeat' to 'All'
	for i in LTR/Ty1_copia LTR/Ty3_gypsy LTR\\t LINE Class_II satellite rDNA All\\t smallClusters singlets; do
		echo -en "${i}\t" >> ${file}_sumpropgen.txt
		grep -P "$i" ${file}_sumprop2.txt | awk 'BEGIN{FS="\t"; sum1=0; sum2=0; sum3=0} {sum1+=$2; sum2+=$3; sum3+=$4} END{print sum1 "\t" sum2 "\t" sum3}' >> ${file}_sumpropgen.txt
	done
	rm ${file}_sumprop2.txt
	sed -i 's/\t\t/\t/' ${file}_sumpropgen.txt #two TABs to a single
	sed -i 's/Class_II/DNAtransposons/' ${file}_sumpropgen.txt
	sed -i 's/All/UnclassifiedRepeats/' ${file}_sumpropgen.txt
	#Add header
	echo -e "Repeat\tReads\tProportion\tMb" > header.txt
	cat header.txt ${file}_sumprop.txt > tmp && mv tmp ${file}_sumprop.txt
	cat header.txt ${file}_sumpropsimple.txt > tmp && mv tmp ${file}_sumpropsimple.txt
	cat header.txt ${file}_sumpropgen.txt > tmp && mv tmp ${file}_sumpropgen.txt
	rm header.txt
	#Read stats
	echo -e "All reads\t${all}" >> ${file}_stat.txt
	echo -e "Top clusters reads\t${topclust}" >> ${file}_stat.txt
	echo -e "Minor clusters reads\t${minorclust}" >> ${file}_stat.txt
	echo -e "Plastid reads\t${plastid}" >> ${file}_stat.txt
	echo -e "Mitochondrial reads\t${mitochondrial}" >> ${file}_stat.txt
	echo -e "Nuclear reads\t${allnuclear}" >> ${file}_stat.txt
	echo -e "Top clusters nuclear reads\t${topclustnuc}" >> ${file}_stat.txt
	echo -e "Singlets\t${singlets}" >> ${file}_stat.txt
	echo -e "Genome size (entered)\t${gs} Mb" >> ${file}_stat.txt
	echo "Finished for ${sample}..."
	echo
	cd ..
done

#Summary table of all samples
echo -e "\nSummarizing all samples..."
#create headers
echo -e "Creating headers..."
echo -e "Species" > RE_stat.txt
echo -e "All reads" >> RE_stat.txt
echo -e "Top clusters reads" >> RE_stat.txt
echo -e "Minor clusters reads" >> RE_stat.txt
echo -e "Plastid reads" >> RE_stat.txt
echo -e "Mitochondrial reads" >> RE_stat.txt
echo -e "Nuclear reads" >> RE_stat.txt
echo -e "Top clusters nuclear reads" >> RE_stat.txt
echo -e "Singlets" >> RE_stat.txt
echo -e "Genome size (entered)" >> RE_stat.txt
echo -e "Repeat\tMb" > RE_sumpropsimpleMb.txt
echo -e "Repeat" > RE_sumpropgenMb.txt
echo -e "LTR/Ty1_copia" >> RE_sumpropgenMb.txt
echo -e "LTR/Ty3_gypsy" >> RE_sumpropgenMb.txt
echo -e "LTR" >> RE_sumpropgenMb.txt
echo -e "LINE" >> RE_sumpropgenMb.txt
echo -e "DNAtransposons" >> RE_sumpropgenMb.txt
echo -e "satellite" >> RE_sumpropgenMb.txt
echo -e "rDNA" >> RE_sumpropgenMb.txt
echo -e "UnclassifiedRepeats" >> RE_sumpropgenMb.txt
echo -e "smallClusters" >> RE_sumpropgenMb.txt
echo -e "singlets" >> RE_sumpropgenMb.txt
#loop over all samples and append results
echo -e "Merging data..."
for file in $(cat listGS.txt | cut -f1); do
	echo -e "...${file}"
	#1. Stat
	echo $file > ${file}_stat.txt
	cat ${file}/${file}_stat.txt | cut -f2 >> ${file}_stat.txt
	paste RE_stat.txt ${file}_stat.txt > tmp && mv tmp RE_stat.txt
	rm ${file}_stat.txt
	#2. Sumpropsimple
	#get only 1st (names) and 4th (Mb) columns, change Mb to species name
	cat ${file}/${file}_sumpropsimple.txt | cut -f1,4 | sed "s/Mb/${file}/" > ${file}_sumpropsimpleMb.txt
	#join the files, add 'NA' if no match, TAB as separator
	join -a1 -a2 -e "NA" -o auto -t $'\t' <(head -1 RE_sumpropsimpleMb.txt; sort <(sed -n '2,$p' RE_sumpropsimpleMb.txt)) <(head -1 ${file}_sumpropsimpleMb.txt; sort <(sed -n '2,$p' ${file}_sumpropsimpleMb.txt)) > tmp && mv tmp RE_sumpropsimpleMb.txt
	rm ${file}_sumpropsimpleMb.txt
	#3. Sumpropgen
	cat ${file}/${file}_sumpropgen.txt | cut -f4 | sed "s/Mb/${file}/" > ${file}_sumpropgenMb.txt
	paste RE_sumpropgenMb.txt ${file}_sumpropgenMb.txt > tmp && mv tmp RE_sumpropgenMb.txt
	rm ${file}_sumpropgenMb.txt
done
#Remove 2nd column from RE_sumpropsimpleMb.txt (containg initial NAs only)
cat RE_sumpropsimpleMb.txt | cut --complement -f2 > tmp && mv tmp RE_sumpropsimpleMb.txt
#Rename 'All' to 'UnclassifiedRepeats' in RE_sumpropsimpleMb.txt
sed -i 's/All/UnclassifiedRepeats/' RE_sumpropsimpleMb.txt

#Make a plot using R
echo -e "\nCreating plots..."
module add r/4.4.0-gcc-10.2.1-ssuwpvb
R --slave -f REplot.R

echo -e "\nFinished RE results summary in ${folder}...\n"
