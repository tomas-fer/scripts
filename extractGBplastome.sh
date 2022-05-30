#EXTRACT sequences (in FASTA format) from full plastome provided in GenBenbank format
#- all features (CDS, tRNA, rRNA), separated to exons
#- all sequences among them (introns, spacers)

#PREPARE reference for HybPhyloMaker
#- modify names
#- delete regions shorter than $short (defined below!)

#USAGE: ./extractGBplastome.sh plastomeFile.gb
#------------------------------------------------------------------------------------
#Tomas Fer, 2018, 2020
#Department of Botany, Charles University, Prague, Czech Republic
#tomas.fer@natur.cuni.cz
#------------------------------------------------------------------------------------

#KNOWN BUGS:
#1. limited working in case the '/gene' descriptor is missing in CDS description (like, e.g., in MK460223)
#   in this case this script requires 'product.txt' for proper modification to gene name ('/gene' descriptor is added)
#   however, some '/product' descriptions matching appropriate gene names have to be added/modified to 'product.txt'
#   known problem: 'RNA polymerase beta subunit' or 'RNA polymerase beta' describes rpoB, rpoC1 and rpoC2 - manual check necessary!!!
#   known problem: some CDS have '/product' description as 'no product string in file'
#   also, if tRNAs and rRNAs are not annotated in GB file, some of the non-coding regions extracted by this script include them
#2. Does not properly work for plastomes of parazitic plants, e.g. HG514460
#   in this case only remaining functional genes are extracted (those annotated as CDS)
#   non-coding regions then contain non-coding parts as well as non-functional genes (i.e. genes with '/pseudo')
#   tRNAs are extracted both functional and pseudo


#regions shorter than this value will be removed
short=200

name=`echo $1 | cut -d'.' -f1`
#get just the sequence
echo -e "\nExtracting sequence..."
#take everything after 'ORIGIN', remove numbers, remove spaces, remove '//', remove CRs, remove LFs
sed -e '1,/ORIGIN/d' $1 | sed 's/[0-9]*//g' | sed 's/ //g' | sed 's/\/\///' | sed 's/\x0D$//' | tr -d '\n' > ${name}.fasta

#Modify GB file (in order to put join feature information on a single line)
#delete everything between ',' and any number [0-9], i.e. putting the join description of some CDS on a single line
perl -0777 -pe 's/(,)\s*([0-9])/$1$2/g' $1 > $1.modif #'perl -0777' = do not any splitting

# Modify GB file in case there are no '/gene' descriptors for CDS
if [ "$(grep "/gene" $1 | wc -l)" -lt 20 ]; then
	echo "Converting '/products' to '/genes'..."
	# modify long names of '/product' (the case of 'ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit')
	perl -0777 -pe 's/(oxygenase)\s*(large)/$1 $2/g' $1.modif > tmp && mv tmp $1.modif
	# change '/product' to gene name according to product.txt
	cat product.txt | while read -r a b; do
		c=$(echo "$b" | sed 's/\//\\\//g') #change '/' to '\/' before sed on the next line
		sed -i "s/\"$c\"/\"$a\"/" $1.modif
	done
	# delete EOL on lines with 'note'
	sed -i '/note/ { N; s/\n// }' $1.modif
	# remove unnecessary lines
	sed -i '/\/note/ d' $1.modif
	sed -i '/\/trans_splicing/ d' $1.modif
	sed -i '/\/codon_start/ d' $1.modif
	sed -i '/\/transl/ d' $1.modif
	sed -i '/locus_tag/ d' $1.modif
	sed -i '/\/exception/ d' $1.modif
	# change '/product' to '/gene'
	sed -i 's/\/product/\/gene/' $1.modif
fi

echo "Extracting features..."
#select lines containing CDS, tRNA and one following line (the line with gene name)
#CDS is surrounded by spaces ( CDS ) to avoid capturing the string from AA translation part
grep -A1 -E " CDS |tRNA " $1.modif > resultCDS.txt
#add '--'' as last line (to mimic grep output after later file merging)
#select lines containing rRNA and two following lines (the line with gene name)
echo "--" >> resultCDS.txt
#remove spaces in gene names if any
sed -i '/gene/s/ //g' resultCDS.txt
#replace two or more spaces by ',', remove lines starting with ',/locus' or ',/gene' or '/old_locus_tag' and select lines containing rRNA and one following line
#(this is to remove other descriptors and get product' descriptor as the next line to rRNA)
sed 's/ \{2,\}/,/g' $1.modif | sed '/^,\/locus/ d' | sed '/^,\/gene/ d' | sed '/^,\/old_locus_tag/ d' | grep -A1 -E ",rRNA" --no-group-separator | sed "s/,rRNA/--\n,rRNA/" | sed 1d > resultR.txt
rm $1.modif
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
#delete last character on last line (if it is an excessive ',')
sed -i '$s/,$//' result.txt
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

#modifying reference for the use with HybPhyloMaker
echo "Modifying reference for HybPhyloMaker..."
#1. separate to types (tRNA, rRNA, CDS, non, CDS+non)
grep "tRNA$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_tRNA.fasta
grep "rRNA$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_rRNA.fasta
grep "CDS$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_CDS.fasta
grep "non$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_non.fasta
grep -e "CDS$" -e "non$" -A1 --no-group-separator ${name}_reference.fasta > ${name}_reference_CDSnon.fasta
#2. modify names
#first sed command replaces last occurrence of '_' by 'zzz'
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_tRNA.fasta > ${name}_reference_tRNAHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_rRNA.fasta > ${name}_reference_rRNAHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_CDS.fasta > ${name}_reference_CDSHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_non.fasta > ${name}_reference_nonHPM.fasta
sed 's/\(.*\)_/\1zzz/; s/_/exon/3; s/_/exon/3; s/-/XXX/; s/\./yy/' ${name}_reference_CDSnon.fasta > ${name}_reference_CDSnonHPM.fasta
#3. remove regions shorter than 200 bp
grep '^[^>].\{200\}' -B1 --no-group-separator ${name}_reference_tRNAHPM.fasta > ${name}_reference_tRNAHPM200.fasta
grep '^[^>].\{200\}' -B1 --no-group-separator ${name}_reference_rRNAHPM.fasta > ${name}_reference_rRNAHPM200.fasta
grep '^[^>].\{200\}' -B1 --no-group-separator ${name}_reference_CDSHPM.fasta > ${name}_reference_CDSHPM200.fasta
grep '^[^>].\{200\}' -B1 --no-group-separator ${name}_reference_nonHPM.fasta > ${name}_reference_nonHPM200.fasta
grep '^[^>].\{200\}' -B1 --no-group-separator ${name}_reference_CDSnonHPM.fasta > ${name}_reference_CDSnonHPM200.fasta

#statistics
echo "Calculating nr of regions..."
#full reference
echo -e "CDS\nnon\ntRNA\nrRNA" > ${name}_statHeader.txt
grep "CDS$" ${name}_reference.fasta | wc -l > ${name}_statD.txt
grep "non$" ${name}_reference.fasta | wc -l >> ${name}_statD.txt
grep "tRNA$" ${name}_reference.fasta | wc -l >> ${name}_statD.txt
grep "rRNA$" ${name}_reference.fasta | wc -l >> ${name}_statD.txt
paste ${name}_statHeader.txt ${name}_statD.txt > ${name}_statT.txt
rm ${name}_statHeader.txt ${name}_statD.txt

#after short removal
grep "CDS$" ${name}_reference_CDSHPM200.fasta | wc -l > ${name}_statD.txt
grep "non$" ${name}_reference_nonHPM200.fasta | wc -l >> ${name}_statD.txt
grep "tRNA$" ${name}_reference_tRNAHPM200.fasta | wc -l >> ${name}_statD.txt
grep "rRNA$" ${name}_reference_rRNAHPM200.fasta | wc -l >> ${name}_statD.txt
paste ${name}_statT.txt ${name}_statD.txt > ${name}_statF.txt
rm ${name}_statT.txt ${name}_statD.txt
echo -e "type\tfull\tlonger200" > header.txt
cat header.txt ${name}_statF.txt > ${name}_stat.txt
rm header.txt ${name}_statF.txt

echo -e "FINISHED\n"
