#EXTRACT sequences (in FASTA format) from full plastome provided in GenBenbank format
#- all features (CDS, tRNA, rRNA), separated to exons
#- all sequences among them (introns, spacers)

#USAGE: ./extractGBplastome.sh plastomeFile.gb
#------------------------------------------------------------------------------------
#Tomas Fer, 2018
#Department of Botany, Charles University, Prague, Czech Republic
#tomas.fer@natur.cuni.cz
#------------------------------------------------------------------------------------

#KNOWN BUGS:
#- problem with parsing GenBank files if CDS line is on multiple lines

name=`echo $1 | cut -d'.' -f1`
#get just the sequence
echo "Extracting sequence..."
#take everything after 'ORIGIN', remove numbers, remove spaces, remove '//', remove CRs, remove LFs
sed -e '1,/ORIGIN/d' $1 | sed 's/[0-9]*//g' | sed 's/ //g' | sed 's/\/\///' | sed 's/\x0D$//' | tr -d '\n' > ${name}.fasta

echo "Extracting features..."
#select lines containing CDS, tRNA and one following line (the line with gene name)
grep -A1 -E "CDS|tRNA " $1 > resultCDS.txt
#add '--'' as last line (to mimic grep output after later file merging)
#select lines containing rRNA and two following lines (the line with gene name)
echo "--" >> resultCDS.txt
#replace two or more spaces by ',', remove lines starting with ',/locus' or ',/gene' and select lines containing rRNA and one following line
#(this is to remove other descriptors and get product' descriptor as the next line to rRNA)
sed 's/ \{2,\}/,/g' $1 | sed '/^,\/locus/ d' | sed '/^,\/gene/ d' | grep -A1 -E ",rRNA" > resultR.txt
#grep -A2 -E "rRNA " $1 > resultR.txt
sed -i 's/ ribosomal RNA//' resultR.txt
#replace spaces by ','
sed -i 's/ \{1,\}/,/g' resultCDS.txt
sed -i 's/ \{1,\}/,/g' resultR.txt
#delete first ',' at each line
sed -i 's/,//' resultCDS.txt
sed -i 's/,//' resultR.txt
#delete first '/' at each line
sed -i 's/\///' resultCDS.txt
sed -i 's/\///' resultR.txt
#select correct name by deleting other lines
sed -i '/^locus/ d' resultR.txt
#merge both files
cat resultCDS.txt resultR.txt > result.txt
rm resultCDS.txt resultR.txt
#replace EOLs by ','
cat result.txt | tr "\n" "," > tmp && mv tmp result.txt
#delete "(" and ")"
sed -i 's/(//g' result.txt
sed -i 's/)//g' result.txt
#delete 'complement'
sed -i 's/complement//g' result.txt
#add ',' after 'join'
sed -i 's/join/join,/g' result.txt
#replace ',--,' by EOL
sed -i "s/,--,/\n/g" result.txt
#replace '..' by ','
sed -i 's/\.\./,/g' result.txt
#delete 'gene='
sed -i 's/gene=//' result.txt
#delete 'product='
sed -i 's/product=//' result.txt
#remove double quotes
sed -i 's/\"//g' result.txt
#delete last character on last line (which is an excessive ',')
sed -i '$ s/.$//' result.txt
#add newline at the end of the file
sed -i '$a\' result.txt
#replace '-' by 'x'
sed -i 's/-/x/' result.txt

echo "Treating joined features..."
#extract lines containing join
grep join result.txt > join.txt
#calculates number of regions to join
awk -F "," '{ print ( NF - 3 ) / 2 }' join.txt > nrRegions.txt
#change ',' to TABs
cat join.txt | tr "," "\t" > tmp && mv tmp join.txt
#merge two files
paste nrRegions.txt join.txt > nrRegionsJoin.txt

#loop over genes with more parts and separates parts to one line each
cat nrRegionsJoin.txt | while read a b c d e f g h i j; do
	if [ "$a" -eq 2 ]; then
		echo -e "$b,$d,$e,${h}_1" >> joinResult.txt
		echo -e "$b,$f,$g,${h}_2" >> joinResult.txt
	elif [ "$a" -eq 3 ]; then
		echo -e "$b,$d,$e,${j}_1" >> joinResult.txt
		echo -e "$b,$f,$g,${j}_2" >> joinResult.txt
		echo -e "$b,$h,$i,${j}_3" >> joinResult.txt
	fi
done
rm join.txt nrRegions.txt nrRegionsJoin.txt

#remove lines with 'join' from results.txt
sed -i '/join/d' result.txt
#combine result.txt and joinResult.txt
cat result.txt joinResult.txt > finalResult.txt
rm result.txt joinResult.txt
#change ',' to TABs
cat finalResult.txt | tr "," "\t" > tmp && mv tmp finalResult.txt
#sort according 2nd field (i.e., starting position)
sort -n -k2 finalResult.txt > finalResultSorted.txt
rm finalResult.txt

#add intergenic regions
echo "Adding intergenic regions..."
#generate file starting positions of intergenic regions
awk '{ print $3 + 1 "\t" $4}' finalResultSorted.txt > startIntergenic.txt
#get name of the last gene
last=$(tail -1 finalResultSorted.txt | cut -f4)
#add as first line
sed -i "1i 1\t$last" startIntergenic.txt
#merge sorted file and 
paste finalResultSorted.txt startIntergenic.txt > finalPlusIntergenic.txt
rm startIntergenic.txt
#write noncoding regions
awk '{ print "non\t" $5 "\t" $2 - 1 "\t" $6 "-" $4}' finalPlusIntergenic.txt > intergenic.txt
rm finalPlusIntergenic.txt
#get the length of the plastome
length=$(head -n1 $1 | sed 's/ \{1,\}/,/g' | cut -d',' -f3)
#get the end of last coding gene and increase by one (start position of last non-coding region)
lastcodpos=$(tail -n1 finalResultSorted.txt | cut -f3)
lastcodpos=$((lastcodpos+1))
#delete last line
head -n -1 intergenic.txt > tmp && mv tmp intergenic.txt
#create last line
echo -e "non\t$lastcodpos\t$length\t$last-end" > lastline.txt
cat intergenic.txt lastline.txt > tmp && mv tmp intergenic.txt
rm lastline.txt

#merge with coding and non-coding and sort
cat finalResultSorted.txt intergenic.txt | sort -n -k2 > tmp && mv tmp finalResultSorted.txt
rm intergenic.txt

#treat negative-length non-coding regions (established by the script when genes are overlapping)
echo "Treating negative-lengths regions..."
#calculate region lengths and add it as forth column
awk '{ print $1 "\t" $2 "\t" $3 "\t" $3 - $2 + 1 "\t" $4 }' finalResultSorted.txt > tmp && mv tmp finalResultSorted.txt
#print only lines with non-negative region length (4th column > 0)
awk '{ if ($4 > 0) print $0 }' finalResultSorted.txt > tmp && mv tmp finalResultSorted.txt

#create reference
echo "Creating reference..."
x=1
cat finalResultSorted.txt | while read a b c d e; do
	echo -e ">${x}_${x}_${e}_${a}" >> ${name}_reference.fasta
	cut -c $b-$c ${name}.fasta >> ${name}_reference.fasta
	x=$((x+1))
done
echo -e "feature\tstart\tstop\tlength\tname" > header.txt
cat header.txt finalResultSorted.txt > ${name}_plastomeFeaturePositions.txt
rm header.txt finalResultSorted.txt
