#Download FASTQ files from BaseSpace using API
#All fastq.gz files from a project are downloaded
#Access token must be in the text file (token_header.txt) whitch contain one line text:                                            #
#header = "x-access-token: <your-token-here>
#----------------------------------------------------------------------------
#Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2019
#tomas.fer@natur.cuni.cz
#----------------------------------------------------------------------------

#Specify the project ID here
projectID=70305241

#Samples from projects
#get info about samples
echo -e "\nGetting info about samples in the project ${projectID}..."
curl -L -J --config ./token_header.txt https://api.basespace.illumina.com/v1pre3/projects/${projectID}/samples?Limit=1000 2>/dev/null > JSONproject.txt
#get sample numbers
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
exit
cat filesList.txt | parallel 'curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/{}/content -O'

#Check file sizes of downloaded files (if they match sizes stated by BaseSpace)
#Download incorrectly downloaded files
echo -e "\nChecking whether file sizes are correct..."
cat fileTable.txt | sed '1d' | while read line; do
	fileName=$(cut -f1 <<< $line) #file name
	fileSizeDown=$(stat -c %s `cut -f1 <<< $line`) #file size of the downloaded file
	fileSizeBS=$(cut -f2 <<< $line) #file size extracted from BaseSpace JSON
	fileBS=$(cut -f3 <<< $line) #file BaseSpace ID
	if [[ ${fileSizeDown} -eq ${fileSizeBS} ]]; then
		echo ${fileName} size OK
	else
		echo -e "${fileName} size incorrect. Downloading again...\n"
		mv ${fileName} ${fileName}.bak 2>/dev/null #make a backup of the wrongly downloaded file
		until [[ $(stat -c %s `cut -f1 <<< $line` 2>/dev/null) -eq ${fileSizeBS} ]]; do
			curl -L -J --config token_header.txt https://api.basespace.illumina.com/v1pre3/files/${fileBS}/content -O
			echo
		done
	fi
done

echo -e "\nFinished doownloading" `cat filesList.txt | wc -l` "files from the project '$expName' (ID: $projectID)"
exit

#Download information about specific run
#Specify the run ID here
runID=76397324

#Run details (total yield, number of clusters, total reads, total reads PF...)
echo "Getting info about run ${runID}"
#Get info about the run
curl -L -J --config ./token_header.txt https://api.basespace.illumina.com/v1pre3/runs/${runID} 2>/dev/null > JSONrun.txt

grep -Po '"YieldTotal":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > yield.txt
grep -Po '"PercentPf":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > percent.txt
grep -Po '"ReadsTotal":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > reads.txt
grep -Po '"ReadsPfTotal":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > readsPF.txt
grep -Po '"Clusters":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > clusters.txt
grep -Po '"ClustersPf":.*?[^\\]"' JSONrun.txt | awk -F\" '{print $3}' | sed 's/[:,]//g' > clustersPF.txt

#make a table
echo -e "TotalYield\tClusters\tClustersPF\tReadsTotal\tReadsTotalPF\tPercentPF" > runTable.txt
paste yield.txt clusters.txt clustersPF.txt reads.txt readsPF.txt percent.txt >> runTable.txt
rm yield.txt clusters.txt clustersPF.txt reads.txt readsPF.txt percent.txt

echo -e "Summary of the run $runID is in 'runTable.txt'"
