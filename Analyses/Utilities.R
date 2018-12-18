
# Running non-parametric comparisons of values 
#@param variable_stack: input vector of the variables to be compared "stacked" on top of each other
#@param number_of_variables: number of variables

VarDist.test<-function(variable_stack, number_of_ratios){
  # stack the three ratios underneath each other, add a factor for the ANOVA
  
  
  ratios1=as.data.frame(cbind(variable_stack, rep(1:number_of_ratios,each=length(variable_stack)/number_of_ratios)))
  colnames(ratios1)=c("Variables","Factor")
  as.factor(ratios1$Factor)
  
  #run nonparametric comparison of means
  kruskal_test<-kruskal.test(Variables~as.factor(Factor), data=ratios1 )
  
  #only if this is significant, consider pairwise wilcoxon test with bonferroni holmes adjustment for multiple comparisons; this is a paired test so it explicitly compares data points from within each species
  pairwise_wilcox<-pairwise.wilcox.test(ratios1$Variables, ratios1$Factor, paired=TRUE, p.adjust.method = "BH")
  
  return(list(kruskal_test,pairwise_wilcox))
  
}









#Determining the best model under different phylogenetic correlation structures and original vs. Grafen-transformed branch lengths

#@param formula: formula of interest (e.g. Tlength~Wt_g), , data are the data   (assuming 
#@param tree: phylogenetic tree, matched to tree in the data preparation step
#@param taildata: data, matched to tree in the data preparation step










#output is output$AIC, then use the name for the highest weighted model to find output$model_summaries$bestmodel
  
AIC=list()
summary_names=c("corPagel","corBrownian", "corGrafen","corPagel_GrafScale","corBrownian_GrafScale", "corGrafen_GrafScale" )

find.best.model<-function(formula, Tree, taildata) {
  
    Pagel=summary (gls(formula, correlation=corPagel(1,phy=Tree), data=taildata))
    Brownian=summary (gls(formula, correlation=corBrownian(1,phy=Tree), data=taildata))
    Grafen=summary (gls(formula, correlation=corGrafen(1,phy=Tree), data=taildata))
    Pagel_GrafScale=summary (gls(formula, correlation=corPagel(1,phy=compute.brlen(Tree)), data=taildata))
    Brownian_GrafScale=summary (gls(formula, correlation=corBrownian(1,phy=compute.brlen(Tree)), data=taildata))
  
    Grafen_GrafScale=summary (gls(formula, correlation=corGrafen(1,phy=compute.brlen(Tree)), data=taildata))
  
    #list models, give names  
    models=list(Pagel, Brownian,Grafen, Pagel_GrafScale,Brownian_GrafScale,Grafen_GrafScale)
    names(models)<-as.vector(summary_names)
  
    #determine AICs
    AIC=list()
  
    for (i in 1:length(names(models))){
    AIC[[i]]=models[[i]]$AIC
    }
  
  #name AICs with correct correlation structures
  names(AIC)<-as.vector(summary_names)
  
  #determine the most likely model
  AICs=unlist(AIC)
  AICmin=AICs-min(AICs)
  W=exp(-0.5*AICmin)/sum(exp(-0.5*AICmin))
  
  #prepare output of most likely model and searchable relevant summary
  all_info=(list(W,models))
  names(all_info)=c("AIC_Weights", "Model_summaries")
  
  #"return" is required to get the  out of the function
  return(all_info)
  
}

 
 

