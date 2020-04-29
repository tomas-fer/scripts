# scripts
various smaller scripts for phylogenomics etc
---
__astralScoring.r__
  * processing scored tree from ASTRAL (with the option -t 4)
  * creates a PNG image of the rooted ASTRAL species tree with pie charts on branches
  * requirement: treeio

__quartetsampling.r__
  * processing trees from quartet sampling (Pease et al., 2018)
  * output trees are (1) modified with sed&grep and (2) plotted in R to resemble trees in publications
  * requirements: devtools, phyloch, phytools

__monophyly.R__
  * testing whether multiple groups are monophyletic in multiple trees
  * loop over is.monophyletic function
  * requirements: ape, phytools

__cpDNA_mapping.sh__
  * BWA mapping of filtered reads from HybPhyloMaker to cpDNA reference
  * consensus call with kindel
  * combine sequences to single FASTA file
  * requirements: bwa, samtools, kindel v.0.1.4

__extractGBplastome.sh__
  * extract sequences (in FASTA format) from full plastome provided in GenBank format: (1) all features (CDS, tRNA, rRNA), separated to exons, (2) all sequences among them (introns, spacers)
  
__BaseSpaceDownload.sh__
  * download FASTQ files from Illumina BaseSpace using API
  * reports basic information about Illumina run
  * requirements: GNU parallel, curl
