####Panphlan Part: Strains based on gene repertoire####
jaccard_gene=as.data.frame(as.matrix(vegdist(t(gene)),index="jaccard"))
jaccard_gene$cohort=NA
for(i in 1:(nrow(jaccard_gene))){
  jaccard_gene[grep("X",row.names(jaccard_gene)),255]=c("HIV")
  jaccard_gene[grep("LL",row.names(jaccard_gene)),255]=c("DAG")
  jaccard_gene[grep("HV",row.names(jaccard_gene)),255]=c("500FG")
}
cluster_ja=hclust(as.dist(jaccard_gene[,1:254]), method = "complete", members = NULL)
par(cex.axis=2,cex=0.5,mar = c(15,4,2,2),bg = 'white')
plot(hBB_ja)

library(dendextend)

types <- factor(jaccard_gene$cohort,levels = c("HIV","DAG","500FG"))
n_types <- length(unique(jaccard_gene$cohort))
cols <- colorspace::rainbow_hcl(n_types, c = 70, l  = 50)

colors_to_use <- cols[types]
hBB_ja=as.dendrogram(hclust(as.dist(jaccard_gene[,1:254]), method = "complete", members = NULL),hang = -1, labels_cex = 0.5)
labels_colors(hBB_ja) <- colors_to_use[order.dendrogram(hBB_ja)]
par(cex.axis=2,cex=0.5,mar = c(5,0,0,0),bg = 'white')
library(circlize)
ID_col_fun = colorRamp2(c(1, 2,3), c("#2a9d8f", "#48cae4","#e76f51"))
color=ID_col_fun(as.numeric(as.factor(jaccard_gene$cohort[cluster_ja$order])))

hBB_ja %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1) %>%  # node point size
  set("leaves_col",color) %>% # node point color
  plot

####draw heatmap of gene repertoire####
gene8=t(gene)
gene8=merge(gene8,cut_ja, 
            by.x = "row.names",
            by.y = "Row.names",
            all=F)
gene8$Strain=as.factor(gene8$Strain)
gene8$Strain=as.numeric(gene8$Strain)
gene8$ID=as.factor(gene8$ID)
gene8$ID=as.numeric(gene8$ID)
row.names(gene8)=gene8$Row.names
gene8$Row.names=NULL
gene8=gene8[row.names(jaccard_gene),]

library(ComplexHeatmap)
library(circlize)
col_fun<-colorRamp2(breaks=c(0,1),
                    colors=c("#ffccd5","#a4133c"))
strain_cluster_col_fun = colorRamp2(c(1, 2, 3), c("#116f51", "#d4008f", "#b50000")) 
ID_col_fun = colorRamp2(c(1, 2, 3), c("#48cae4", "#2a9d8f", "#e76f51"))#500FG:3 DAG:2 HIV:1

top = HeatmapAnnotation(
  Strain_cluster = anno_simple(gene8$Strain, col = strain_cluster_col_fun),
  Cohort=anno_simple(gene8$ID, col = ID_col_fun),
  annotation_name_side = "right")

row_anno=merge(as.data.frame(t(gene9))[,1,drop=F],DA_gene_2clusters,
               by.x = "row.names",
               by.y = "gene",
               all=T)
row_anno=na.omit(row_anno)
row.names(row_anno)=row_anno$Row.names
row_anno$Row.names=NULL
row_anno$V1=NULL
col_fun2<-colorRamp2(breaks=c(min(row_anno[,2]),1e-100,0.05),
                     colors=c("red","black","#ffffff"))
row_ha=HeatmapAnnotation(`P value`=row_anno[,2],
                         col = list(`P value` = col_fun2),which="row")
HH=Heatmap(t(gene8[ ,1:4260]), name = "gene", col = col_fun,
           column_names_gp = gpar(fontsize =10),
           row_names_gp = gpar(fontsize = 1),
           column_names_rot=45,
           cluster_rows=T,
           cluster_columns =as.dendrogram(hclust(as.dist(jaccard_gene[,1:254]), method = "complete", members = NULL),hang = -1, labels_cex = 0.5) , 
           top_annotation=top)
HH+row_ha

####enrichment####
cut_ja=as.data.frame(cutree(cluster_ja,k=3))
cut_ja=merge(cut_ja,ID4, 
             by.x = "row.names",
             by.y = "row.names",
             all=F)
cut_ja$group2=NA
cut_ja$group2[which(cut_ja$ID==c("500FG")|cut_ja$ID==c("DAG"))]=c("healthy")
cut_ja$group2[which(cut_ja$ID==c("HIV"))]=c("HIV")
colnames(cut_ja)[2]=c("Strain")
fisher.test(table(cut_ja[-which(cut_ja$Strain==3),c(2,11)]))$p.value
table(cut_ja[,c(6,7)])
fisher.test(table(cut_ja[,c(6,7)]))$p.value
#draw
ggplot(cut_ja, 
       aes(x = Strain, fill=ID)) + 
  geom_bar(position = "stack")

####GO enrichment in 2 clusters####
library(clusterProfiler)# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html#ora-algorithm
library(enrichplot)
library(ggplot2)
####GO enrichment analysis in HIV-related strain####
colnames(DA_gene_2clusters0.05)[1]=c("NR90")
ego <- enricher(gene= DA_gene_2clusters0.05$NR90,# differential gene families
                TERM2GENE= annotation4[c("GO","NR90")], # all gene families
                pAdjustMethod= "BH",
                pvalueCutoff= 0.05,
                qvalueCutoff= 0.1)
yy=ego@result

####GO enrichment analysis in control-related strain####
colnames(DA_gene_2clusters0.052)[1]=c("NR90")
ego2 <- enricher(gene= DA_gene_2clusters0.052$NR90,# differential gene
                 TERM2GENE= annotation4[c("GO","NR90")], # all gene families
                 pAdjustMethod= "BH",
                 pvalueCutoff= 0.05,
                 qvalueCutoff= 0.05)
yy2=ego2@result
####draw####
ego2=pairwise_termsim(ego)

emapplot(ego2, cex_label_category=.8, cex_line=.5) + 
  scale_fill_continuous(low = "#e06663", high = "#327eba", name = "p.adjust",
                        guide = guide_colorbar(reverse = TRUE, order=1), trans='log10') 

cnetplot(ego,showCategory = 5,foldChange=DA_gene_2clusters0.05$beta)

####DA gene between clusters####
gene9=merge(gene8[,1:4260],cut_ja,
            by.x="row.names",
            by.y="Row.names",
            all=F)
gene9=gene9%>% `rownames<-`(.[,1]) %>% select(-Row.names)
logistic_regression_gene=function(data){
  results=data.frame(gene=NA,
                     cluster.p.value=NA,
                     cluster.BH_FDR=NA,
                     beta=NA)
  data$Strain=as.factor(data$Strain)
  data$GEN=as.factor(data$GEN)
  for(i in grep("Uni",colnames(data))){
    model1=glm(data[,i]~GEN+AGE+log10_read+Strain,data = data)
    results[i,1]=colnames(data)[i]
    sum1=summary(model1)
    results[i,2]=sum1$coefficients[5,4] ## Cluster2 vs.Sub-cluster1_p.value
    results[i,4]=sum1$coefficients[5,1] ## Cluster2 vs.Sub-cluster1_p.value
  }
  results[,3]= p.adjust(results[,'cluster.p.value'],method = "BH")
  results<-results
}
DA_gene_2clusters=logistic_regression_gene(gene9[!gene9$Strain==3,])
DA_gene_2clusters0.05=DA_gene_2clusters[DA_gene_2clusters$cluster.BH_FDR<0.05&DA_gene_2clusters$beta>0,]
DA_gene_2clusters0.052=DA_gene_2clusters[DA_gene_2clusters$cluster.BH_FDR<0.05&DA_gene_2clusters$beta<0,]



####IL-10####
#residules of cytokine production corrected for gender+age+read counts+SO+Num-P
strain_IL10=merge(anti[,5,drop=F],cut_ja[!cut_ja$Strain==3,],
               by.x = "row.names",
               by.y = "Row.names",
               all=F)
strain_IL10$group2=as.factor(strain_IL10$group2)
strain_IL10=merge(strain_IL10,t(HIV_noantibio[grep("s__Prevotella_copri",row.names(HIV_noantibio)),,drop=F]),
               by.x = "Row.names",
               by.y = "row.names",
               all=F)
strain_IL10$s__Prevotella_copri=as.numeric(strain_IL10$s__Prevotella_copri)
strain_IL10$Strain=as.factor(strain_IL10$Strain)
strain_IL10=merge(strain_IL10,pheno[,c(19,20)],
                  by.x = "Row.names",
                  by.y = "row.names",
                  all=F)

####sex enrichment####
fisher.test(table(strain_IL10[,c(3,13)]))$p.value
table(strain_IL10[,c(3,13)])
fisher.test(matrix(c(4,7,5,68),nrow = 2,ncol = 2))
####IL-10####
data1=strain_IL10[strain_IL10$Strain==1,]
row.names(data1)=data1$Row.names
residule1=as.data.frame(lapply(data1[,2,drop=F],function(x){resid(lm(x~AGE+GEN+log10_read+SO+`Num-P`,data=data1))}))

data2=strain_IL10[strain_IL10$Strain==2,]
row.names(data2)=data2$Row.names
residule2=as.data.frame(lapply(data2[,2,drop=F],function(x){resid(lm(x~AGE+log10_read+SO+`Num-P`,data=data2))}))

data1=merge(data1,residule1,
            by.x = "row.names",
            by.y = "row.names",
            all=F)
data2=merge(data2,residule2,
            by.x = "row.names",
            by.y = "row.names",
            all=F)
IL10_result1=lm(Pam3Cys_PBMC_24h_IL.10.x~s__Prevotella_copri,data=data1)
IL10_result2=lm(Pam3Cys_PBMC_24h_IL.10.y~s__Prevotella_copri,data=data2)

sum1=summary(IL10_result1)
sum2=summary(IL10_result2)

#heterogeneity test
studlab=c("1","2")
TE=c(sum1$coefficients[2,1],sum2$coefficients[2,1])
seTE=c(sum1$coefficients[2,2],sum2$coefficients[2,2])
inDf <- data.frame(studlab, TE, seTE)
library(meta)
meta_stain_1and2 <- metagen(TE = TE,
                            seTE = seTE,
                            studlab = studlab,
                            data = inDf,
                            title = "1 vs.2")
p.heterogeneity=meta_stain_1and2$pval.Q

alldata=rbind(data1,data2)
alldata$Row.names=NULL

p=ggplot(data = alldata, mapping = aes(x =s__Prevotella_copri, y =`Pam3Cys_PBMC_24h_IL.10.y`))+
  geom_point(aes(color=Strain))+
  geom_smooth(data =alldata, mapping = aes(x =s__Prevotella_copri, y = `Pam3Cys_PBMC_24h_IL.10.y`,fill=Strain,color=Strain),method = "lm", formula = y ~ x)
library(ggExtra)
ggMarginal(p, type = "violin",groupFill=T )
ggMarginal(p, type = "density",groupFill=T )

####demographics for 2 strains####
#continues varialbe
demographics_strain=merge(cut_ja[!cut_ja$Strain==3,1:2],pheno2_inv,
                          by.x = "Row.names",by.y = "row.names",
                          all=F)
demographics_strain$Strain=factor(demographics_strain$Strain,levels=c("2","1"))
wilcox=function(data){
  results=data.frame(pheno=NA,
                     wilcoxon.p.value=NA,
                     wilcoxon.BH_FDR=NA)
  for(i in c(3:12,14:19,22)){
    model=wilcox.test(data[,i]~Strain,data = data)
    results[i,1]=colnames(data)[i]
    results[i,2]=model$p.value
  }
  results[,3]= p.adjust(results[,'wilcoxon.p.value'],method = "BH")
  results<-results
}
demographics_strain_result=wilcox(demographics_strain)
#Categorical varialbe
fisher.test(table(demographics_strain[,c(2,13)]))
fisher.test(table(na.omit(demographics_strain[!demographics_strain$SO==1,c(2,21)])))
fisher.test(table(na.omit(demographics_strain[,c(2,23)])))
####draw
p1=ggplot (demographics_strain, aes(Strain,`Num-P`,color=Strain))+
  geom_violin(alpha = 0.3,width=0.6)+
  geom_jitter(size=3,width = 0.08, height = 0.1,alpha=0.3)+
  xlab(label = "Strain") + ylab(label = "Num-P") +
  theme(legend.text = element_text(size = 8, face = "bold"),legend.position = c(0.08,0.9))+
  scale_x_discrete(breaks= c("1","2"),
                   labels = c("Health-related","HIV-related"))

p2=ggplot (demographics_strain, aes(Strain,`CD4-recovery-abs`,color=Strain))+
  geom_violin(alpha = 0.3,width=0.6)+
  geom_jitter(size=3,width = 0.08, height = 0.1,alpha=0.3)+
  xlab(label = "Strain") + ylab(label = "Num-P") +
  theme(legend.text = element_text(size = 8, face = "bold"),legend.position = c(0.08,0.9))+
  scale_x_discrete(breaks= c("1","2"),
                   labels = c("Health-related","HIV-related"))

p3=ggplot (demographics_strain, aes(Strain,`CD4 counts`,color=Strain))+
  geom_violin(alpha = 0.3,width=0.6)+
  geom_jitter(size=3,width = 0.08, height = 0.1,alpha=0.3)+
  xlab(label = "Strain") + ylab(label = "Num-P") +
  theme(legend.text = element_text(size = 8, face = "bold"),legend.position = c(0.08,0.9))+
  scale_x_discrete(breaks= c("1","2"),
                   labels = c("Health-related","HIV-related"))

p7=ggplot (demographics_strain, aes(Strain,`CA-DNA`,color=Strain))+
  geom_violin(alpha = 0.3,width=0.6)+
  geom_jitter(size=3,width = 0.08, height = 0.1,alpha=0.3)+
  xlab(label = "Strain") + ylab(label = "CA-DNA")+
  scale_x_discrete(breaks= c("1","2"),
                   labels = c("Health-related","HIV-related"))
p8=ggplot (demographics_strain, aes(Strain,`CA-RNA`,color=Strain))+
  geom_violin(alpha = 0.3,width=0.6)+
  geom_jitter(size=3,width = 0.08, height = 0.1,alpha=0.3)+
  xlab(label = "Strain") + ylab(label = "CA-RNA")+
  scale_x_discrete(breaks= c("1","2"),
                   labels = c("Health-related","HIV-related"))

p9=ggplot (demographics_strain, aes(Strain,`CA-RNA/CA-DNA`,color=Strain))+
  geom_violin(alpha = 0.3,width=0.6)+
  geom_jitter(size=3,width = 0.08, height = 0.1,alpha=0.3)+
  xlab(label = "Strain") + ylab(label = "CA-RNA")+
  scale_x_discrete(breaks= c("1","2"),
                   labels = c("Health-related","HIV-related"))

p7+p8+p9

library(patchwork)
(p1+p2+p3)/(p4+p5+p6)

demographics_strain$GEN=as.factor(demographics_strain$GEN)
demographics_strain$SO=as.factor(demographics_strain$SO)
demographics_strain$RAI=as.factor(demographics_strain$RAI)

p4=ggplot(demographics_strain, 
       aes(x = Strain, fill=GEN)) + 
  geom_bar(position = "stack")

p5=ggplot(demographics_strain, 
       aes(x = Strain, fill=SO)) + 
  geom_bar(position = "stack")

p6=ggplot(demographics_strain, 
       aes(x = Strain, fill=RAI)) + 
  geom_bar(position = "stack")

####uniref90 annotation####
Glutamate_5_kinase_uniref90<-read.table("~/Downloads/Glutamate_5_kinase.tsv",sep="\t",quote = "",row.names = NULL,header = T,check.names = F,fill = T)




