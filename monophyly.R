#TESTING MONOPHYLY OF PREDEFINED GROUPS IN GENE TREES
#multiple loop over is.monophyletic()
#by Tomas Fer, 2016
#-------------------------------------------------------------

#load libraries
library (phytools)
library (ape)

#DEFINE MONOPHYLETIC GROUPS OF INTEREST
#Groups for exons Zingiberales_subset5_structured_MORE_GROUPS
Musineae=c("Musella-lasiocarpa_S132","Heliconia-aurantiaca_S345","Heliconia-pastazae_S349","Strelitzia-nicolai_S351","Ensete-superbum_S375","Heliconia-bourgaeana_S347","Phenakospermum-guyanense_S142","Musa-acuminata_cds","Heliconia-lingulata_S348","Ravenala-madagascariensis_S181","Orchidantha-chinensis_S352","Musa-balbisiana_cds","Heliconia-bihai_S346","Heliconia-waoraniana_S350")
allExclMusaceae=c("Phrynium-imbricatum_S133","Calathea-crotalifera_S401","Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Heliconia-aurantiaca_S345","Heliconia-pastazae_S349","Strelitzia-nicolai_S351","Sarcophrynium-sweinfurthianum_S353","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Heliconia-bourgaeana_S347","Phenakospermum-guyanense_S142","Canna-liliiflora_S402","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Heliconia-lingulata_S348","Ravenala-madagascariensis_S181","Orchidantha-chinensis_S352","Stachyphrynium-sinense_S377","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125","Heliconia-bihai_S346","Heliconia-waoraniana_S350")
Zing=c("Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125")
Cost=c("Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185")
Mar=c("Phrynium-imbricatum_S133","Calathea-crotalifera_S401","Stachyphrynium-sinense_S377","Sarcophrynium-sweinfurthianum_S353")
Strel=c("Strelitzia-nicolai_S351","Phenakospermum-guyanense_S142","Ravenala-madagascariensis_S181")
Helic=c("Heliconia-aurantiaca_S345","Heliconia-pastazae_S349","Heliconia-bourgaeana_S347","Heliconia-lingulata_S348","Heliconia-bihai_S346","Heliconia-waoraniana_S350")
Mus=c("Musella-lasiocarpa_S132","Ensete-superbum_S375","Musa-acuminata_cds","Musa-balbisiana_cds")
CostZing=c("Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125")
CanMar=c("Canna-liliiflora_S402","Phrynium-imbricatum_S133","Calathea-crotalifera_S401","Stachyphrynium-sinense_S377","Sarcophrynium-sweinfurthianum_S353")
HelicCanMarCostZing=c("Heliconia-aurantiaca_S345","Heliconia-pastazae_S349","Heliconia-bourgaeana_S347","Heliconia-lingulata_S348","Heliconia-bihai_S346","Heliconia-waoraniana_S350","Canna-liliiflora_S402","Phrynium-imbricatum_S133","Calathea-crotalifera_S401","Stachyphrynium-sinense_S377","Sarcophrynium-sweinfurthianum_S353","Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125")
LowStrel=c("Strelitzia-nicolai_S351","Phenakospermum-guyanense_S142","Ravenala-madagascariensis_S181","Orchidantha-chinensis_S352")
LowStrelHelic=c("Strelitzia-nicolai_S351","Phenakospermum-guyanense_S142","Ravenala-madagascariensis_S181","Orchidantha-chinensis_S352","Heliconia-aurantiaca_S345","Heliconia-pastazae_S349","Heliconia-bourgaeana_S347","Heliconia-lingulata_S348","Heliconia-bihai_S346","Heliconia-waoraniana_S350")
LowStrelMus=c("Strelitzia-nicolai_S351","Phenakospermum-guyanense_S142","Ravenala-madagascariensis_S181","Orchidantha-chinensis_S352","Musella-lasiocarpa_S132","Ensete-superbum_S375","Musa-acuminata_cds","Musa-balbisiana_cds")
MusHelic=c("Musella-lasiocarpa_S132","Heliconia-aurantiaca_S345","Heliconia-pastazae_S349","Ensete-superbum_S375","Heliconia-bourgaeana_S347","Musa-acuminata_cds","Heliconia-lingulata_S348","Musa-balbisiana_cds","Heliconia-bihai_S346","Heliconia-waoraniana_S350")
CostZingCanMar=c("Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125","Canna-liliiflora_S402","Sarcophrynium-sweinfurthianum_S353","Phrynium-imbricatum_S133","Stachyphrynium-sinense_S377","Calathea-crotalifera_S401")
ZingCan=c("Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125","Canna-liliiflora_S402")
CostZingCan=c("Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125","Canna-liliiflora_S402")
CostCanMar=c("Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Canna-liliiflora_S402","Sarcophrynium-sweinfurthianum_S353","Phrynium-imbricatum_S133","Stachyphrynium-sinense_S377","Calathea-crotalifera_S401")
allExclHeliconiaceae=c("Musella-lasiocarpa_S132","Ensete-superbum_S375","Musa-acuminata_cds","Musa-balbisiana_cds","Phrynium-imbricatum_S133","Calathea-crotalifera_S401","Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Strelitzia-nicolai_S351","Sarcophrynium-sweinfurthianum_S353","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Phenakospermum-guyanense_S142","Canna-liliiflora_S402","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Ravenala-madagascariensis_S181","Orchidantha-chinensis_S352","Stachyphrynium-sinense_S377","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125")
allExclLowStrel=c("Heliconia-aurantiaca_S345","Heliconia-pastazae_S349","Heliconia-bourgaeana_S347","Heliconia-lingulata_S348","Heliconia-bihai_S346","Heliconia-waoraniana_S350","Musella-lasiocarpa_S132","Ensete-superbum_S375","Musa-acuminata_cds","Musa-balbisiana_cds","Phrynium-imbricatum_S133","Calathea-crotalifera_S401","Tapeinochilos-ananassae_S182","Costus-spectabilis_S269","Costus-curvibracteatus_S246","Tamijia-flagellaris_S22","Hedychium-ellipticum_S245","Renealmia-polypus_S66","Hornstedtia-bella_S144","Sarcophrynium-sweinfurthianum_S353","Chamaecostus-cuspidatus_S191","Costus-talbotii_S189","Costus-scaber_S305","Globba-marantina_S33","Zingiber-citriodorum_S171","Alpinia-nigra_S81","Canna-liliiflora_S402","Paracostus-paradoxus_S266","Monocostus-uniflorus_S187","Costus-phaeotrichus_S306","Siphonochilus-aethiopicus_S130","Curcuma-flaviflora_S315","Pleuranthodium-racemigerum_S135","Amomum-tsaoko_S136","Stachyphrynium-sinense_S377","Costus-tonkinensis_S268","Dimerocostus-argenteus_S190","Costus-afer_S185","Siphonochilus-decorus_S24","Monolophus-appendiculatus_S289","Siliqamomum-tonkinense_S56","Alpinia-oblongifolia_S125")

#Musaceae-related groups only (comment if you want to test all groups mention above)
groupsToTest=list(Musineae, allExclMusaceae, Zing, Cost, Mar, Strel, Helic, Mus, CostZing, CanMar, HelicCanMarCostZing, LowStrel, LowStrelHelic, LowStrelMus, MusHelic, CostZingCanMar, ZingCan, CostZingCan, CostCanMar, allExclHeliconiaceae, allExclLowStrel)
groupsToTestNames=c("Musineae","allExclMusaceae","Zing","Cost","Mar","Strel","Helic","Mus","CostZing","CanMar","HelicCanMarCostZing","LowStrel","LowStrelHelic","LowStrelMus","MusHelic","CostZingCanMar","ZingCan","CostZingCan","CostCanMar","allExclHeliconiaceae","allExclLowStrel")

#define trees (all files with *.tre in current directory) - this command will store all tree names in 'trees_files'
trees_files <- dir(pattern="*.result")

#function reading newick tree and evaluating monophyly of species in a group (tips = "group")
monophyletic <- function(file, monolist) {
  print (file)
  #read newick tree
  tree = read.newick(file)
  #put all species to 'alltips'
  alltips <- tree$tip.label
  #in 'comparelist' will be only those species from 'monolist' that are present in a tree
  comparelist <- alltips[alltips %in% monolist]
  #if nr. of species in comparelist is at least two function is.monophyletic is called
  if (length(comparelist) > 1) {
    mono <- is.monophyletic(phy = tree, tips = comparelist)
    #otherwise function returns "NA"
  } else {
    mono <- NA
  }
  return(c(mono))
}

#put names of trees to ismonophyl
ismonophyl <- trees_files
#make dataframe (trees names in the first column)
ismonophyl <- data.frame(tree=trees_files)

x <- 0
#loop over all groups
for(i in groupsToTest) {
  x <- x + 1
  #apply monophyly test of a group to all trees
  ismonophyl[[x+1]] <- lapply(trees_files, monophyletic, monolist = i)
}

#add column names to the matrix
colnames(ismonophyl) <- c("tree", groupsToTestNames)
#make matrix
ismonophyl <- as.matrix(ismonophyl)
write.csv(ismonophyl, file="monophylyResults.txt", quote=FALSE, row.names=FALSE)
#***SCRIPT END***
