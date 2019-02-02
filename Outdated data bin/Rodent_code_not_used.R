#Rodent tree code in case it's of interest
#The rodent tree needs some processing from the file supplied by the publication 
rodTree_raw=read.nexus("../Data/Smissen_Rowe_2018_Rodent_Tree.nex")


rod=rod[-which(is.na(rod$Tl_mm)),]
rod=rod[-which(is.na(rod$Locomotor_use)),]

#remove catalogue numbers
rodTree=rodTree_raw

split_tips=list()
new_tips=list()
tips_for_processing=rodTree$tip.label; 

for (i in 1:length(tips_for_processing)){
  
  split_tips[i]=strsplit(tips_for_processing[i], "_")
  
  new_tips[i]=paste (split_tips[[i]][1],split_tips[[i]][2], sep="_" )
}

newtips=unlist(new_tips)

#replace old tip labels with new ones
rodTree$tip.label=newtips

#remove duplicate species names from tree
rodTree=drop.tip(rodTree, which(duplicated(rodTree$tip.label)))

#find typo positions, fix
rodTree$tip.label[grep("Notomys_mithcellii", rodTree$tip.label)]="Notomys_mitchellii"
rodTree$tip.label[grep("Conilurus_pencillatus", rodTree$tip.label)]=c("Conilurus_penicillatus")


#Rodents: 
namecheck=name.check(rodTree, rod)
rodTree_synch=drop.tip(rodTree,namecheck$tree_not_data)
#in rodents there are also a few species that are not in the phylogeny, unfortunately:
rod_synch=rod[-which(rownames(rod) %in% namecheck$data_not_tree), ] 
name.check(rodTree_synch,rod_synch)

#Murines
murine_data=rod[which(rod$Order =="Murine"),]
murine_tree=drop.tip(rodTree_synch, name.check(rodTree_synch,murine_data)$tree_not_data) 

#not sufficient resolution for Rattus



#in analyses file ot make table of two rows of data quality
as.data.frame(rbind(sapply(data_raw[,2:28], function(x) sum((!is.na(x))/nrow(data_raw))*100),sapply(rodraw[,2:28], function(x) sum((!is.na(x)))/nrow(rodraw))*100))
rownames(DataQuality)=c("Marsupials", "Rodents")


