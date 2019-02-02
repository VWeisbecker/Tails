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



#Rodents

##return the min/max values to their un-logged form
ratio_Tl=c(10^(rodFull$Tl_max)/10^(rodFull$Tl_min))
ratio_Bl=c(10^(rodFull$Bl_max)/10^(rodFull$Bl_min))
ratio_Wt=c(10^(rodFull$Wt_max)/10^(rodFull$Wt_min))

##Visualizing means of the ratios 
boxplot(ratio_Tl, ratio_Bl,ratio_Wt)

ratiosRod=c(ratio_Tl, ratio_Bl, ratio_Wt)

##The VarDist.test code is in utilities. It runs a Kruskal-Wallis test and then a pairwise wilcoxon test on the ratios
MinMaxtest(ratiosRod,3)

boxplot(ratio_Tl, ratio_Bl,ratio_Wt)

##Is there a difference between Rattus and the old endemics.

ratio_Tl_Rat<-ratio_Tl[ - (grep("Rattus_leucopus", rodFull$Species):grep("Rattus_fuscipes", rodFull$Species))]
ratio_Bl_Rat<-ratio_Bl[ - (grep("Rattus_leucopus", rodFull$Species):grep("Rattus_fuscipes", rodFull$Species))]
ratio_Wt_Rat<-ratio_Wt[ - (grep("Rattus_leucopus", rodFull$Species):grep("Rattus_fuscipes", rodFull$Species))]

ratiosRat=c(ratio_Tl_Band, ratio_Bl_Band, ratio_Wt_Band)
MinMaxtest(ratiosRat,3)

par(mfrow=c(1,2))
boxplot(ratio_Tl, ratio_Bl,ratio_Wt, xlab="ratios with Rattus")
boxplot(ratio_Tl_Rat, ratio_Bl_Rat,ratio_Wt_Rat, xlab="ratios without Rattus")







#Rodents:

RodMales<-cbind(rodFull$Tl_Av_M,rodFull$Bl_Av_M,rodFull$Wt_Av_M, rep(1,length(rodFull$Tl_Av_M)))
RodFemales<-cbind(rodFull$Tl_Av_F,rodFull$Bl_Av_F,rodFull$Wt_Av_F, rep(2,length(rodFull$Tl_Av_F)))
RodAverages<-cbind(rodFull$Tl_mm, rodFull$Bl_mm,rodFull$Wt_g, rep(3,length(rodFull$Bl_mm)))

RodCombine=rbind(RodMales, RodFemales,RodAverages)
RodCombine=as.data.frame(RodCombine)
colnames(RodCombine) <-c("Tl","Bl","Wt","Sex")
RodCombine$Sex<-as.factor(RodCombine$Sex)


##Are there differences in slope of tail length and body lenght/weight according to locomotor use, sex, or averages between sexes?
###check for appropriate distribution of model residuals
plot(lm(RodCombine$Tl~RodCombine$Bl*RodCombine$Sex))

summary(aov(lm(RodCombine$Tl~RodCombine$Bl*RodCombine$Sex)))

###No interaction (i.e. slopes are the same), therefore dropping interaction
summary(aov(lm(RodCombine$Tl~RodCombine$Bl+RodCombine$Sex)))

##Same with body weight
plot(lm(RodCombine$Tl~RodCombine$Wt*RodCombine$Sex))
summary(aov(lm(RodCombine$Tl~RodCombine$Wt*RodCombine$Sex)))
##No interaction (i.e. slopes are the same), therefore dropping interaction
summary(aov(lm(RodCombine$Tl~RodCombine$Wt+RodCombine$Sex)))




#Rodents with Rattus

rodTBW=find.best.model(Tl_mm~Bl_mm*Wt_g,rodTree,rod);rodTW$AIC_Weights
#Checking appropriate distribution of model residuals
qqnorm (rodTBW$Models$Brownian)
anova(rodTBW$Models$Brownian)

###dropping the interaction 
rodTBWNoInter=find.best.model(Tl_mm~Bl_mm+Wt_g,rodTree,rod);rodTW$AIC_Weights
#Checking appropriate distribution of model residuals; this one is not ideal
qqnorm (rodTBWNoInter$Models$Brownian_GrafScale)
#significant but with small effect size
anova(rodTBWNoInter$Models$Brownian_GrafScale, type="marginal")

#Rodents without Rattus

MurTBW=find.best.model(Tl_mm~Bl_mm*Wt_g,murine_tree,murine_data);MurTW$AIC_Weights
#Checking appropriate distribution of model residuals
qqnorm (MurTBW$Models$Brownian_GrafScale)
anova(MurTBW$Models$Brownian_GrafScale, type="marginal")

###dropping the interaction 
MurTBWNoInter=find.best.model(Tl_mm~Bl_mm+Wt_g,murine_tree,murine_data);MurTW$AIC_Weights
#Checking appropriate distribution of model residuals; this one is not ideal
qqnorm (MurTBWNoInter$Models$Brownian_GrafScale)
#significant but with small effect size; but higher level of significance without Rattus
anova(MurTBWNoInter$Models$Brownian_GrafScale, type="marginal")


