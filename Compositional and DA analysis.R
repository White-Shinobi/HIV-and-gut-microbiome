####data####

#relative abundance
#HIV:HIV1,HH
#DMP:DAG1,DD
#HIV and DMP: HIV_DAG3,beta_pathway

####packages####
library(vegan)
library(igraph)
library(circlize)

####species level####
####alpha diversity
HIV_shannon=as.data.frame(diversity(t(HIV1),index="shannon"))
DAG3_shannon=as.data.frame(diversity(t(DAG1),index="shannon"))
wilcox_test_shannon=wilcox.test(as.numeric(HIV_shannon[,1]),as.numeric(DAG3_shannon[,1]),exact = T, correct = FALSE)
####beta diversity
bray=as.matrix(vegdist(HIV_DAG3),index="bray")
PERMANOVA=adonis(bray~BMI+log10_read+ID,ID9[row.names(HIV_DAG3),],permutations = 1000)#put covariate before the factors you're interested in!!!

####pathway level####
####alpha diversity####
HIV_shannon2=as.data.frame(diversity(t(HH),index="shannon"))
DAG3_shannon2=as.data.frame(diversity(t(DD),index="shannon"))
wilcox_test_shannon2=wilcox.test(as.numeric(HIV_shannon2[,1]),as.numeric(DAG3_shannon2[,1]),exact = T, correct = FALSE)
####beta diversity####
beta_pathway_bray=as.matrix(vegdist(t(beta_pathway)),index="bray")
PERMANOVA_pathway=adonis(beta_pathway_bray~BMI+log10_read+ID,ID9[row.names(beta_pathway_bray),],permutations = 1000)

####Differential abundance analysis between HIV and DMP####

#clr_transformation
do_clr_externalWeighting = function(interest_matrix, core_matrix){
  #######how to deal with zero
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  #geometric mean=exp(mean(log(x)))
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 2, gm_mean)
  
  #do transformation
  data_prepared = rbind(Gmean_core,interest_matrix)
  data_transformed = apply(data_prepared,2,function(x){
    log(x / x[1])[-1]
  })
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}
####invrank####
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

####linear regression: find out DA spp####
linear_regression=function(data,abundace_DAG,abundace_HIV){
  results=data.frame(microbe=NA,variables=NA,HIVCohort_p.value=NA,HIVCohort_BH_FDR=NA,HIVCohort_estimate=NA,
                     Cohor_CI2.5=NA,Cohor_CI97.5=NA,
                     mean_DMP=NA,mean_HIV=NA)
  for(i in 1:(ncol(data)-5)){
    dag=abundace_DAG[row.names(abundace_DAG)%in%colnames(data)[i],]
    hiv=abundace_HIV[row.names(abundace_HIV)%in%colnames(data)[i],]
    data$ID=as.factor(data$ID)
    data$smoke=as.factor(data$smoke)
    model=lm(data[,i]~ID+log10_read+BMI+smoke,data = data)
    sum1=summary(model)
    CI=confint(model, "IDHIV",level=0.95) ##cal CI of ID:HIV
    results[i,1]=colnames(data)[i]
    results[i,2]=c("Cohort+log10_read+BMI+smoke")
    results[i,3]=sum1$coefficients[2,4] ## Cohort_p.value
    results[i,5]=sum1$coefficients[2,1] ## Cohort_estimate
    results[i,6]=CI[1,1] ## Cohor_CI2.5
    results[i,7]=CI[1,2] ## Cohor_CI97.5
    results[i,8]=mean(as.numeric(dag[1,]),na.rm = T)
    results[i,9]=mean(as.numeric(hiv[1,]),na.rm = T)
  }
  results[,4]= p.adjust(results[,'HIVCohort_p.value'],method = "BH")
  results<-results
}
#input data:clr first, then filter 20% and invrank transformed
DA_result=linear_regression(all5,DAG1,HIV1)
DA_result0.05=DA_result[DA_result$HIVCohort_BH_FDR<0.05,]

write.table(DA_result,"1.tsv",sep="\t",row.names = T,col.names = T)
sum(DA_result0.05$HIVCohort_estimate>0)#57 species enriched in PLHIV
sum(DA_result0.05$HIVCohort_estimate<0)#19 species depleted in PLHIV

DApathway_result=linear_regression(all4_inv,DD,HH)
DApathway_result0.05=DApathway_result[DApathway_result$HIVCohort_BH_FDR<0.05,]
DApathway_result0.05$`log2(HIV/DMP)`=log2(DApathway_result0.05$mean_HIV/DApathway_result0.05$mean_DMP)
DApathway_result0.05=merge(DApathway_result0.05,annotation3,
                           by.x = "microbe",
                           by.y = "Bacterial pathway",
                           all=F)
sum(DApathway_result0.05$HIVCohort_estimate>0)#87 species enriched in PLHIV
sum(DApathway_result0.05$HIVCohort_estimate<0)#76 species enriched in PLHIV

DAclass_result=linear_regression(all_class_final,DD,HH)
DAclass_result0.05=DAclass_result[DAclass_result$HIVCohort_BH_FDR<0.05,]
sum(DAclass_result0.05$HIVCohort_estimate>0)#13 species enriched in PLHIV
sum(DAclass_result0.05$HIVCohort_estimate<0)#10 species enriched in PLHIV

####draw heatmap of DA species####
library(tidyverse)
heat_DA=all5%>%select(DA_result0.05$microbe[order(DA_result0.05$HIVCohort_p.value,decreasing = F)][1:30],ID)
colnames(heat_DA)=sapply(strsplit(as.character(colnames(heat_DA)),"s\\__"),"[",2)
colnames(heat_DA)=gsub("\\_"," ",as.character(colnames(heat_DA)))
colnames(heat_DA)[31]=c("ID")
heat_DA$ID=as.numeric(heat_DA$ID)
library(circlize)
col_fun<-colorRamp2(breaks=c(min(heat_DA[,-31]),0,max(heat_DA[,-31])),
                    colors=c("#03045e","white","#ef476f"))
ID_col_fun = colorRamp2(c(1, 2), c("#4057AB","#E26971")) #DAG:1 HIV:2

library(ComplexHeatmap)
top = HeatmapAnnotation(
  Cohort=anno_simple(heat_DA$ID, col = ID_col_fun),
  annotation_name_side = "right")

Heatmap(t(heat_DA[,-31]),
        row_dend_side = "right",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        name = "Color key", 
        col = col_fun,
        column_names_gp = gpar(fontsize =10),
        row_names_gp = gpar(fontsize = 10),
        column_names_rot=45,
        top_annotation=top)

####draw network of DA pathways####

#from--to
graph_path=DApathway_result0.05[,c(1,3,5,13)]
colnames(graph_path)[c(1,4)]=c("from","to")
graph_path=graph_path[,c(1,4,2,3)]
#metadata for all nodes：all side nodes and the central nodes
#Central nodes
center_node=as.data.frame(unique(graph_path$to))
center_node=merge(center_node,DAclass_result0.05[,c(1,5)],
                  by.x = "unique(graph_path$to)",
                  by.y = "microbe",
                  all=T)
colnames(center_node)[1]=c("vertices")
center_node$label=center_node$vertices
center_node$shape=c("square")
center_node$color=NA
center_node$color[center_node$HIVCohort_estimate<0]=c("#3d52a0")
center_node$color[center_node$HIVCohort_estimate>0]=c("#d8676e")
center_node$color[is.na(center_node$HIVCohort_estimate)]=c("white")
center_node$HIVCohort_p.value=10^(-5)
colnames(center_node)[2]=colnames(beside_node)[2]
#Side nodes
beside_node=graph_path[,c(1,3,4)]
beside_node$label=c("")
colnames(beside_node)[1]=c("vertices")
beside_node$shape=c("circle")
#Annotate nodes colors

col_fun<-colorRamp2(breaks=c(min(beside_node$HIVCohort_estimate),0,max(beside_node$HIVCohort_estimate)),
                    colors=c("#3d52a0","white","#d8676e"))   #023047","white","#ffb703"
beside_node$color=col_fun(beside_node$HIVCohort_estimate)
beside_node=beside_node[,colnames(center_node)]
all_node_meta=rbind(center_node,beside_node)
#Draw the figure
g <- graph_from_data_frame(d = graph_path,vertices=all_node_meta,directed=F)
#set the size for all nodes
size=-log10(V(g)$`HIVCohort_p.value`)*0.7
size[!is.finite(size)]=mean(size[is.finite(size)])+1
V(g)$size <- size
#set the label for all nodes
V(g)$label=V(g)$label
#set the shape for all nodes
V(g)$shape=V(g)$shape
#set the color for all nodes
V(g)$color=V(g)$color
#plot
plot(g,layout = layout_with_fr,vertex.label.cex=1,vertex.label.font=3) #pdf:9;9






