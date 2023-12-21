#Calculate number of reads/proportion of all identified repeat types in RepeatExplorer cluster annotation file
#Based on annotations in the 5th (automatic) or 7th (final) column
#Annotation must follow the standards of RepeatExplorer, i.e. always start with 'All' etc.
#Proportion is calculated based on number of nuclear reads (i.e. after cp/mt reads removal)
#Size is corrected to the genome size in Mb/1C (provide an integer, i.e. whole number)

#Tomas Fer, 2022
#tomas.fer@natur.cuni.cz
#Version: 0.3 (2022-07-08)

#Usage: ./makeREsummary.sh annotationtype filename genomeSizeMb1C
#e.g. ./makeREsummary.sh a CLUSTER_TABLE.csv 1361 #summary of automatic annotation (i.e., 5th column)
#e.g. ./makeREsummary.sh f CLUSTER_TABLE.csv 1361 #summary of final annotation (i.e., 7th column)

#Output files:
#filename_{automatic,final}_annotation.txt - nr. reads and annotation of all top clusters (columns 4 and 5/7 from filename)
#filename_stat.txt - read statistics
#filename_sumprop.txt - summary of reads/proportions/sizes(Mb) of all repeat types
#filename_sumpropsimple.txt - the same as the previous but with simplified repeat names
#filename_sumpropgen.txt - summary of reads/proportions/sizes(Mb) of basic repeat types

#Check if input file exists
if [ ! -f "$2" ]; then
	echo "error: input file does not exist" >&2; exit 3
fi
#Check if genome size is an integer
re='^[0-9]+$'
if ! [[ $3 =~ $re ]] ; then
	echo "error: GS is not an integer" >&2; exit 3
fi

#Check if annotation type is 'f' or 'a'
if ! [[ $1 =~ ^f$ || $1 =~ ^a$ ]] ; then
	echo "error: annotation type should be 'f' or 'a'" >&2; exit 3
fi

#Get appropriate annotation column
if [[ $1 =~ ^a$ ]] ; then
	annot=automatic
	#get 4th and 5th columns only, i.e., read numbers and automatic annotation; remove lines with TABs only
	awk '{ print $4"\t"$5"\t" }' $2 | sed '/^[[:blank:]]/d' > annotation.txt
elif [[ $1 =~ ^f$ ]] ; then
	annot=final
	#get 4th and 7th columns only, i.e., read numbers and final annotation; remove lines with TABs only
	awk -F'\t' '{ print $4"\t"$7"\t" }' $2 | sed '/^[[:blank:]]/d' > annotation.txt
fi

#Check if annotations exist
if [ -z "$(awk '{ print $2 }' annotation.txt | grep "All")" ]; then
	echo "error: annotation missing or not in the correct format" >&2; exit 3
fi

#Initial screen output
echo
echo "File: $2"
echo "Annotation: $annot"
echo "Genome size (1C): $3 Mb"
echo

#Remove path (i.e. everything before the last '/') and then the suffix from file name
file=$(sed 's@.*/@@' <<< $2 | cut -d'.' -f1)
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
all=$(grep "Number_of_analyzed_reads" $2 | awk '{ print $2 }')
plastid=$(grep "All/organelle/plastid" ${file}_summary.txt | awk '{ print $2 }')
[ -z "$plastid" ] && plastid=0 #set 'plastid' to '0' if empty
mitochondrial=$(grep "All/organelle/mitochondria" ${file}_summary.txt | awk '{ print $2 }')
[ -z "$mitochondrial" ] && mitochondrial=0 #set 'mitochondrial' to '0' if empty
allnuclear=$(( $all - $plastid - $mitochondrial ))
singlets=$(grep "Number_of_singlets" $2 | awk '{ print $2 }')
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
	mb=$(echo -e "scale=8;$3 * ($perc / 100)" | bc)
	echo -e "${perc}\t${mb}" >> ${file}_summary3.txt
done
mv annotation.txt ${file}_${annot}_annotation.txt
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
echo -e "Genome size (entered)\t${3} Mb" >> ${file}_stat.txt
echo "Finished..."
echo
