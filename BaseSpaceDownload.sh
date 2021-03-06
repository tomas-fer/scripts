#Download FASTQ files from BaseSpace using API
#All fastq.gz files from a project are downloaded
#Optionally, some run statistics is reported
#Access token must be in the text file (token_header.txt) containing one line of text:
#header = "x-access-token: <your-token-here>"
#projectID (and runID) must be specified within the script (see below)
#REQUIRES: GNU parallel, curl
#
#----------------------------------------------------------------------------
#Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2020
#tomas.fer@natur.cuni.cz
#----------------------------------------------------------------------------

#Specify the project ID here
projectID=

#Test if prerequisites are installed
for i in curl parallel; do
	if ! [ -x "$(command -v $i)" ]; then
		echo "Please install '$i' first" && exit
	fi
done

#Test if result files exist
for i in sampleTable.txt fileTable.txt filesList.txt samplesList.txt JSONproject.txt JSONsamples.txt; do
	if [[ -f $i ]]; then
		echo The file '$i' already exists. Delete it or rename before running this script again... && exit
	fi
done

#Samples from projects
#Get info about samples
echo -e "\nGetting info about samples in the project ${projectID}..."
curl -L -J --config ./token_header.txt https://api.basespace.illumina.com/v1pre3/projects/${projectID}/samples?Limit=1000 2>/dev/null > JSONproject.txt
#Check whether the project exists
if grep -q Error "JSONproject.txt"; do
	echo -e "\nYou are not permitted to access the project '$projectID'\n" && exit
fi
#Get sample numbers
grep -Po '"Href":.*?[^\\]",' JSONproject.txt | grep "/samples" | awk -F\" '{print $4}'| awk -F\/ '{print $3}' > samplesList.txt
grep -Po '"SampleId":.*?[^\\]",' JSONproject.txt | awk -F\" '{print $4}' > sampleID.txt
grep -Po '"LibraryName":.*?[^\\]",' JSONproject.txt | awk -F\" '{print $4}' > libName.txt
grep -Po '"TotalReadsPF":.*?[^\\]",' JSONproject.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > readsPF.txt
expName=$(grep -Po '"ExperimentName":.*?[^\\]"' JSONproject.txt | awk -F\" '{print $4}' | head -n1)

echo "There are" `cat samplesList.txt | wc -l` "samples in the project '$expName'"

#Make samples table
echo -e "SampleID\tName\tReadsPF\tBaseSpaceID" > sampleTable.txt
paste sampleID.txt libName.txt readsPF.txt samplesList.txt >> sampleTable.txt
rm sampleID.txt libName.txt readsPF.txt

#Get file numbers
echo -e "\nGetting info about files in the project ${projectID}..."
for i in $(cat samplesList.txt); do
	#download information about files for particular sample
	curl -L -J --config ./token_header.txt https://api.basespace.illumina.com/v1pre3/samples/${i}/files?Extensions=gz 2>/dev/null > JSONsamples.txt
	#get 'Id', display only them and add it to the IDs list
	grep -Po '"Id":.*?[^\\]",' JSONsamples.txt | awk -F\" '{print $4}' >> filesList.txt
	#get sizes
	grep -Po '"Size":.*?[^\\]"' JSONsamples.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> filesSize.txt
	#get file names
	grep -Po '"Path":.*?[^\\]"' JSONsamples.txt | awk -F\" '{print $4}' >> filesNames.txt
done

echo "There are" `cat filesList.txt | wc -l` "files"

#Make table
echo -e "FileName\tSize\tBaseSpaceID" > fileTable.txt
paste filesNames.txt filesSize.txt filesList.txt >> fileTable.txt
rm filesNames.txt filesSize.txt

#Download individual files using parallel
echo -e "\nDownloading fastq.gz files..."

cat filesList.txt | parallel 'curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/{}/content -O'

#Check file sizes of downloaded files (if they match sizes stated by BaseSpace)
#Download incorrectly downloaded files
echo -e "\nChecking whether file sizes are correct..."
cat fileTable.txt | sed '1d' | while read line; do
	fileName=$(awk '{ print $1 }' <<< $line) #file name
	fileSizeDown=$(stat -c %s `awk '{ print $1 }' <<< $line`) #file size of the downloaded file
	fileSizeBS=$(awk '{ print $2 }' <<< $line) #file size extracted from BaseSpace JSON
	fileBS=$(awk '{ print $3 }' <<< $line) #file BaseSpace ID
	if [[ ${fileSizeDown} -eq ${fileSizeBS} ]]; then
		echo ${fileName} size OK
	else
		echo -e "${fileName} size incorrect. Downloading again...\n"
		mv ${fileName} ${fileName}.bak 2>/dev/null #make a backup of the wrongly downloaded file
		until [[ $(stat -c %s `awk '{ print $1 }' <<< $line` 2>/dev/null) -eq ${fileSizeBS} ]]; do
			curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/${fileBS}/content -O
			echo
		done
	fi
done

echo -e "\nFinished downloading" `cat filesList.txt | wc -l` "files from the project '$expName' (ID: $projectID)"
exit #Remove if you want to continue with getting info about a run

#Download information about specific run
#Specify the run ID here
runID=

#Run details (total yield, number of clusters, total reads, total reads PF...)
echo "Getting info about run ${runID}"
#Get info about the run
curl -L -J --config ./token_header.txt https://api.basespace.illumina.com/v1pre3/runs/${runID} 2>/dev/null > JSONrun.txt
grep -Po '"ExperimentName":.*?[^\\]"' JSONrun.txt | head -n1 | awk -F\" '{print $4}' | sed 's/[:,]//g' > rundata.txt
grep -Po '"PlatformName":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $4}' | sed 's/[:,]//g' >> rundata.txt
grep -Po '"YieldTotal":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
grep -Po '"PercentPf":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
grep -Po '"Clusters":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
grep -Po '"ClustersPf":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,}]//g' >> rundata.txt
grep -Po '"PercentGtQ30":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
grep -Po '"PercentGtQ30R1":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
grep -Po '"PercentGtQ30R2":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' >> rundata.txt
#make a table
echo -e "ExperimentName\nPlatformName\nTotalYield[Gbp]\nClusters\nClustersPF\nPercentPF\nPercentGtQ30\nPercentGtQ30R1\nPercentGtQ30R2" > runHeader.txt
paste runHeader.txt rundata.txt > runTable.txt
rm rundata.txt runHeader.txt
echo -e "Summary of the run $runID is in 'runTable.txt'"
