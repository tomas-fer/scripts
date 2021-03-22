#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=2:mem=16gb:scratch_local=4gb
#PBS -j oe
#PBS -N PhyParts_Astral
#PBS -m abe

#Runs PhyParts and PhyParts PieCharts for provided species tree
#Takes all gene trees and only select trees with outgroup
#If name of the species tree is 'ASTRAL', first computes Astral tree

#Tomas Fer, 2021, tomas.fer@natur.cuni.cz

#SETTINGS
treefile=trees_ml_exons.nwk #file with all gene trees
astraltree=ASTRAL #file with a single species tree, if named 'ASTRAL' the species tree will be created with Astral
source=/storage/brno2/home/${LOGNAME}/phyparts #directory on MetaCentrum with the data
astraljar=astral.5.7.7.jar #name of the Astral JAR file in the source directory (also the 'lib' directory has to be there!)
OUTGROUP="Siliquamomum-tonkinense_S56" #outgroup as it appears in trees
nrpptrees= #use only first XXX gene trees for PhyParts (leave empty if all trees should be used)
phypartsbs=0.5 #parameter '-s' for PhyParts (cut-off for bootstrap support in gene trees)
ppcolors="" #colours for PhyParts pie charts (four named colours separated by spaces, the whole expression must be within quotes), see https://en.wikipedia.org/wiki/Web_colors for colour's names, leave empty for original colours
#ppcolors="CornflowerBlue YellowGreen Orange LightSlateGray"

#Add necessary MetaCentrum modules
module add phyparts-0.0.1
module add newick-utils-13042016

cd $SCRATCHDIR

#Copy data
cp $source/$treefile .
if [[ $astraltree =~ "ASTRAL" ]]; then
	#copy ASTRAL
	cp $source/$astraljar .
	cp -r $source/lib .
	#run ASTRAL
	java -jar $astraljar -i $treefile -o Astral.tre 2> Astral.log
	cp Astral.tre $source
	cp Astral.log $source
	astraltree=Astral.tre
else
	#copy species tree
	cp $source/$astraltree .
fi

#Remove '.tre' from species tree name
sptree=$(cut -d'.' -f1 <<< $astraltree)

#Modify Astral tree
#replace everything after ';' by nothing, i.e., delete last space if necessary
#replace ' ' back to '-' and '_'
sed -i 's/;.*/;/' Astral.tre
sed -i 's/ \([^ ]*\) / \1_/g' ${sptree}.tre #replace every second occurrence of ' ' by '_'
sed -i 's/ /-/g' ${sptree}.tre #replace all spaces by '-'

#Reroot Astral tree with $OUTGROUP
nw_reroot -s ${astraltree} $OUTGROUP > tmp && mv tmp ${astraltree}

#Modify labels in gene tree
sed -i 's/XX/-/g' $treefile
sed -i 's/YY/_/g' $treefile

#Get only rooted gene trees (i.e., trees containing outgroup)
nrgenetreesorig=$(wc -l < $treefile)
grep "$OUTGROUP" ${treefile} > tmp
mv tmp ${treefile}

#Check if there are any gene trees left
if [ $(wc -l < ${treefile}) -eq 0 ]; then
	echo -e "\nThere are no gene trees rooted with ${OUTGROUP}. Exiting..."
	exit 3
fi

#Reroot genestrees with $OUTGROUP
nw_reroot -s ${treefile} $OUTGROUP > tmp && mv tmp ${treefile}

#Put single tree per file to directory 'trees'
mkdir trees
split -a 4 -d -l 1 ${treefile} trees/tree_

nrgenetreesrooted=$(ls trees/tree_* | wc -l)

#Subselect only first 'nrpptrees' if necessary (to decrease running time)
if [ -z "$nrpptrees" ]; then #test whether $nrpptrees is empty
	nrpptrees=$nrgenetreesrooted #if empty set to number of rooted gene trees
fi

if [ $nrpptrees -gt $nrgenetreesrooted ]; then #test whether $nrpptrees is higher than nr. of rooted trees
	nrpptrees=$nrgenetreesrooted #if higher set to number of rooted gene trees
fi

if [ $nrpptrees -lt $nrgenetreesrooted ]; then
	echo -e "\nSubselecting $nrpptrees trees..."
	mkdir trees$nrpptrees
	cd trees
	ls tree_* | head -n $nrpptrees | xargs -I{} cp "{}" ../trees$nrpptrees/
	cd ..
	rm -r trees
	mv trees$nrpptrees trees
fi

#Calculate number of gene trees
nrgenetrees=$(ls trees/tree_* | wc -l)

#Set to 0.5 if $phypartsbs is empty
if [ -z "$phypartsbs" ]; then #test whether $phypartsbs is empty
	echo -e "\n'phypartsbs' is not set to any value, using 0.5..."
	phypartsbs=0.5
fi

#Statisctics
echo -e "\nSpecies tree: ${sptree}.tre" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Gene tree file: ${treefile}" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Nr. gene trees: ${nrgenetreesorig}" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Nr. gene trees rooted with '$OUTGROUP': ${nrgenetreesrooted}" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt
echo -e "Using $nrpptrees gene trees" | tee -a phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt

#Make dir for results
mkdir phyparts_${phypartsbs}
mkdir phyparts_${phypartsbs}/trees_res_${nrpptrees}trees

#Run phyparts
echo -e "\nRunning PhyParts with support cutoff ${phypartsbs}..."
# -a what kind of analysis (0 - concon, 1 - fullconcon, 2 - duplications)
# -d directory of trees
# -m  mapping tree (species tree)
# -o prepend output files with this
# -s support cutoff (only keep things with greater support than the one specified)
# -v include verbose output
echo "phyparts -a 1 -d trees -m ${sptree}.tre -o trees_res -s ${phypartsbs} -v"
phyparts -a 1 -d trees -m ${sptree}.tre -o trees_res -s ${phypartsbs} -v > phyparts.log 2>&1

#Run phypartspiecharts
echo -e "\nRunning PhyParts PieCharts..."
echo -e "Using colours: ${ppcolors}"
# --svg_name name of the resulting graphics
# species tree
# prefix of phypart results (i.e., the same as '-o' option in phyparts)
# number of genetrees
if [ -z "$ppcolors" ]; then
	echo "phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees}"
	phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees}
else
	echo "phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors}"
	phypartspiecharts --svg_name phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg ${sptree}.tre trees_res ${nrgenetrees} --colors ${ppcolors}
fi

#Convert SVG to PDF
module add python36-modules-gcc #adds also cairosvg
cairosvg phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg -o phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.pdf

#Copy results back
cp phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.svg phyparts_${phypartsbs}
cp phyparts_${sptree}_BS${phypartsbs}_${nrpptrees}trees.pdf phyparts_${phypartsbs}
cp phypartsinfo_BS${phypartsbs}_${nrpptrees}trees.txt phyparts_${phypartsbs}
cp trees_res.*node* phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
cp trees_res.*hist* phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
cp phyparts.log phyparts_${phypartsbs}/trees_res_${nrpptrees}trees
cp -r phyparts_${phypartsbs} $source

#Clean scratch
if [[ ! $SCRATCHDIR == "" ]]; then
	rm -rf $SCRATCHDIR/*
fi
