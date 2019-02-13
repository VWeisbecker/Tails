
# Running non-parametric comparisons of values;  
#@param variable_stack: input vector of the variables to be compared "stacked" on top of each other
#@param number_of_variables: number of variables

MinMaxtest<-function(variable_stack, number_of_ratios){
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
summary_names=c("Pagel","Brownian", "Grafen","Pagel_GrafScale","Brownian_GrafScale", "Grafen_GrafScale" )

find.best.model<-function(formula, Tree, taildata) {
  
    Pagel= (gls(formula, correlation=corPagel(1,phy=Tree), data=taildata))
    Brownian= (gls(formula, correlation=corBrownian(1,phy=Tree), data=taildata))
    Grafen= (gls(formula, correlation=corGrafen(1,phy=Tree), data=taildata))
    Pagel_GrafScale= (gls(formula, correlation=corPagel(1,phy=compute.brlen(Tree)), data=taildata))
    Brownian_GrafScale= (gls(formula, correlation=corBrownian(1,phy=compute.brlen(Tree)), data=taildata))
    Grafen_GrafScale= (gls(formula, correlation=corGrafen(1,phy=compute.brlen(Tree)), data=taildata))
    
    PagelSum=summary (gls(formula, correlation=corPagel(1,phy=Tree), data=taildata))
    BrownianSum=summary (gls(formula, correlation=corBrownian(1,phy=Tree), data=taildata))
    GrafenSum=summary (gls(formula, correlation=corGrafen(1,phy=Tree), data=taildata))
    Pagel_GrafScaleSum=summary (gls(formula, correlation=corPagel(1,phy=compute.brlen(Tree)), data=taildata))
    Brownian_GrafScaleSum=summary (gls(formula, correlation=corBrownian(1,phy=compute.brlen(Tree)), data=taildata))
    Grafen_GrafScaleSum=summary (gls(formula, correlation=corGrafen(1,phy=compute.brlen(Tree)), data=taildata))
  
    #list model summaries, give names  
    models<-list(Pagel, Brownian,Grafen, Pagel_GrafScale,Brownian_GrafScale,Grafen_GrafScale)
    names(models)<-as.vector(summary_names)
    
    summaries<-list(PagelSum, BrownianSum,GrafenSum, Pagel_GrafScaleSum, Brownian_GrafScaleSum,Grafen_GrafScaleSum)  
    names(summaries)<-as.vector(summary_names)
    #determine AICs
    AIC=list()
  
    for (i in 1:length(names(summaries))){
    AIC[[i]]=summaries[[i]]$AIC
    }
  
  #name AICs with correct correlation structures
  names(AIC)<-as.vector(summary_names)
  
  #determine the most likely model
  AICs=unlist(AIC)
  AICmin=AICs-min(AICs)
  W=exp(-0.5*AICmin)/sum(exp(-0.5*AICmin))
  
  #prepare output of most likely model and searchable relevant summary
  all_info=(list(W,models,summaries))
  names(all_info)=c("AIC_Weights", "Models", "Model_summaries")
  
  #"return" is required to get the  out of the function
  return(all_info)
  
}

 

#A code to predict the y parameters of a given x for pgls (or any other) coefficients; this produces code that can be copied and pasted into the plots to add lines for regressions of specific lengths (avoiding a confusion of lines that extrapolate beyond the actual data points)
#@param model_intercept: a model intercept value, this can be navigated to from a saved model (in pgls, the address is model$coefficients[1])
#@param model_slope: the slope coefficient, in pgls this is model$coefficients[2]
#@param value1, value2: x values for start and end of the line

pgls.line.text <- function (model_intercept, model_slope, value1, value2){
  Y1 <- model_slope*value1+model_intercept
  Y2 <- model_slope*value2+model_intercept
  vector <- c(value1,value2,Y1,Y2)
  
  line_coords <- list(vector)
  
  line <-(paste("lines(x=c(" , line_coords[[1]][1], "," , line_coords[[1]][2],")",",y=c(", line_coords[[1]][3],",", line_coords[[1]][4] , "))", sep="" ))
  
  return(line)
  
}
 

