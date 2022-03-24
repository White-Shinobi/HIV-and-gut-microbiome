####INDEX calculation and comparison among 3 cohorts####
####DI/FI calculation function####
INDEX=function(data){
  up=data[grep("up",data$group),-ncol(data)]
  down=data[grep("down",data$group),-ncol(data)]
  up=up%>%
    mutate_each(type.convert)
  down=down%>%
    mutate_each(type.convert)
  #deal with 0s
  if(any(up==0)) up = up + min(up[up>0])/2
  if(any(down==0)) down = down + min(down[down>0])/2
  #geometric mean of spp enriched/delected in HIV
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_up = apply(up, 2, gm_mean)
  Gmean_down = apply(down, 2, gm_mean)
  Gmean_data=rbind(Gmean_up,Gmean_down)
  #DI and FI=log2(Gmean_up/Gmean_down)
  DI_FI_final = apply(Gmean_data,2,function(x){
    log2(x / x[2])
  })
}

#input:all_index_3cohort,all index data from 3 cohorts
library(FSA)
dunnTest(`Dysbiosis_index`~`ID.x`,all_index_3cohort,method = "bh")
dunnTest(`Function_index`~`ID.x`,all_index_3cohort,method = "bh")
dunnTest(`PB_ratio`~`ID.x`,all_index_3cohort,method = "bh")
wilcox.test(`Dysbiosis_index`~`ID.x`,all_index_3cohort[-which(all_index_3cohort$`ID.x`==c("500FG")),])$p.value
wilcox.test(`PB_ratio`~`ID.x`,all_index_3cohort[-which(all_index_3cohort$`ID.x`==c("500FG")),])
#draw
ggplot(all_index_3cohort,aes(x=all_index_3cohort$`log10(P/B_ratio+min(matrix)/2)`,fill=`ID.x`,alpha = 1/10))+
  geom_density()
#violin
Group=factor(DIall$ID,c("DAG","500FG","HIV"))
ggplot (all_index_3cohort, aes(`ID.x`,`log10(P/B_ratio+min(matrix)/2)`)) + 
  stat_boxplot(geom = "errorbar",width=0.15)+  ###max and min
  geom_boxplot(aes(fill = `ID.x`),width=0.15) +
  geom_violin(alpha = 0.3,width=0.6)+
  theme_classic() +
  xlab(label = "Groups") + ylab(label = "log10(P/B_ratio+min(matrix)/2)") +
  theme(legend.text = element_text(size = 8, face = "bold"))+
  scale_x_discrete(breaks= c("500FG", "DAG","HIV"),
                   labels = c("500FG", "DAG", "HIV"))+
  annotate("text", x = 1.5, y = 2, label = paste("p value=",kruskal.test_DIall$p.value),size=4)+
  theme_linedraw()
#draw
ggplot(all_index_3cohort,aes(x=all_index_3cohort$Function_index,fill=`ID.x`,alpha = 1/10))+
  geom_density()
#violin
Group=factor(DIall$ID,c("DAG","500FG","HIV"))

ggplot (all_index_3cohort, aes(`ID.x`,Function_index)) + 
  stat_boxplot(geom = "errorbar",width=0.15)+  ###max and min
  geom_boxplot(aes(fill = `ID.x`),width=0.15) +
  geom_violin(alpha = 0.3,width=0.6)+
  theme_classic() +
  xlab(label = "Groups") + ylab(label = "FI score") +
  theme(legend.text = element_text(size = 8, face = "bold"))+
  scale_x_discrete(breaks= c("500FG", "DAG","HIV"),
                   labels = c("500FG", "DAG", "HIV"))+
  annotate("text", x = 1.5, y = 2, label = paste("p value=",kruskal.test_DIall$p.value),size=4)+
  theme_linedraw()

####phenotype data####

#phenotype raw data: pheno2

#phenotype corrected for GEN+AGE+log10_read: df  
#phenotype corrected for GEN+AGE+log10_read+SO+Num-P: df2  
#adjusted phenotypes ussing GEN+AGE+log10_read;invrank and then correct
pheno2_inv=pheno2
pheno2_inv[,c(1:9,14:16)]=lapply(pheno2_inv[,c(1:9,14:16)],invrank)
pheno2_inv$SO=as.character(pheno2_inv$SO)
pheno2_inv$SO[which(pheno2_inv$SO==c("HIV_MSM"))]=3
pheno2_inv$SO[which(pheno2_inv$SO==c("HIV_MSW"))]=2
pheno2_inv$SO[which(pheno2_inv$SO==c("HIV_W"))]=1
pheno2_inv$RAI=as.numeric(as.character(pheno2_inv$RAI))

residule_pheno1=lapply(pheno2_inv[c(1:9,14:16)],function(x){resid(lm(x~GEN+AGE+log10_read,pheno2_inv))})
df<-data.frame(matrix(1:143,ncol=1))
for(i in names(residule_pheno1)){
  row<-data.frame(residule_pheno1[[i]])
  colnames(row)=i
  df<-as.data.frame(cbind(df,row))
}
df=df[,-c(1)] 
df=merge(df,pheno2_inv[,c(10,19:21)],
         by.x = "row.names",
         by.y = "row.names",
         all=F)
row.names(df)=df$Row.names
df$Row.names=NULL
df$`HIVRNA load`[df$`HIVRNA load`>0]=1
df$`HIVRNA load`[df$`HIVRNA load`==0]=0
df=df[c(1:12,15,13,14,16)]
for (i in (14:16)){df[,i]=as.numeric(df[,i])}
#additionally adjusted for SO and Num-P
residule_pheno2=lapply(pheno2_inv[c(1:9,14:16)],function(x){resid(lm(x~GEN+AGE+log10_read+SO+`Num-P`,pheno2_inv))})
df2<-data.frame(matrix(1:120,ncol=1))
for(i in names(residule_pheno2)){
  row<-data.frame(residule_pheno2[[i]])
  colnames(row)=i
  df2<-as.data.frame(cbind(df2,row))
}
df2=df2[,-c(1)] 
df2=merge(df2,pheno2[,c(10),drop=F],
          by.x = "row.names",
          by.y = "row.names",
          all=F)
row.names(df2)=df2$Row.names
df2$Row.names=NULL
df2$`HIVRNA load`[df2$`HIVRNA load`>0]=1
df2$`HIVRNA load`[df2$`HIVRNA load`==0]=0

####correlation between phenotypes####
pheno3=pheno2
pheno3$SO=as.character(pheno3$SO)
pheno3$SO[which(pheno3$SO==c("HIV_MSM"))]=3
pheno3$SO[which(pheno3$SO==c("HIV_MSW"))]=2
pheno3$SO[which(pheno3$SO==c("HIV_W"))]=1
pheno3$RAI=as.numeric(as.character(pheno3$RAI))
for (i in names(pheno3)){pheno3[,i]=as.numeric(pheno3[,i])}
library(psych)
cor_13pheno=psych::corr.test(pheno3,method = 'spearman', adjust = 'none')
library(reshape2)
melt_cor_13pheno=melt(cor_13pheno$`r`)
melt_p_13pheno=melt(cor_13pheno$`p`)
melt_pheno=cbind(melt_cor_13pheno,melt_p_13pheno[,3])
colnames(melt_pheno)[3]=c("r")
colnames(melt_pheno)[4]=c("p.value")
melt_pheno$FDR=p.adjust(melt_pheno$p.value,method = "BH")
melt_pheno$final=NA
melt_pheno$final[which(melt_pheno$FDR<0.1)]=melt_pheno$r[which(melt_pheno$FDR<0.1)]
melt_pheno$final[which(is.na(melt_pheno$final))]=0

ggplot(melt_pheno, aes(x=Var1,y=Var2, fill=final))+ geom_tile(aes(fill = final), colour = "transparent")+
  scale_fill_gradient2(low = "#023047",mid="white", high = "#ae2012", 
                       breaks = c(min(melt_pheno$final),0,max(melt_pheno$final)))+
  theme(axis.text.x = element_text(angle =90, vjust=0.5, hjust=1))

####correlation test####
####correlation between shannon index, P/B ration, DI and FI score with HIV-related phenotypes#### 
PB_DI_FI_shannon_phe=merge(all_index_3cohort[,c(1,2,5,6,7)],df,
                           by.x = "Row.names",
                           by.y = "row.names",
                           all=F)
row.names(PB_DI_FI_shannon_phe)=PB_DI_FI_shannon_phe$Row.names
PB_DI_FI_shannon_phe$Row.names=NULL
PB_DI_FI_shannon_phe$`HIVRNA load`=as.numeric(PB_DI_FI_shannon_phe$`HIVRNA load`)
PB_DI_FI_shannon_phe$SO=as.numeric(PB_DI_FI_shannon_phe$SO)
PB_DI_FI_shannon_phe$RAI=as.numeric(PB_DI_FI_shannon_phe$RAI)

correlation_4index=function(x){
  results=data.frame(x=NA,y=NA,
                     p.value=NA,BH_FDR=NA,rho=NA)
  k=0
  for(j in c(5:ncol(x))){
    for(i in c(1:4)){
      model=cor.test(x[,i], x[,j],method = "spearman")
      results[(i+k),1]=colnames(x)[i]
      results[(i+k),2]=colnames(x)[j]
      results[(i+k),3]=model$p.value ## p value second row, 4th col
      results[(i+k),5]=model$estimate ## beta coef of DI 
    }
    k=k+100
  }
  results[,'BH_FDR']= p.adjust(results[,'p.value'],method = "BH")
  results<-na.omit(results)
}
PB_DI_FI_shannon_phe_result1=correlation_4index(PB_DI_FI_shannon_phe)

PB_DI_FI_shannon_phe2=merge(all_index_3cohort[,c(1,2,5,6,7)],df2,
                            by.x = "Row.names",
                            by.y = "row.names",
                            all=F)
row.names(PB_DI_FI_shannon_phe2)=PB_DI_FI_shannon_phe2$Row.names
PB_DI_FI_shannon_phe2$Row.names=NULL
PB_DI_FI_shannon_phe2$`HIVRNA load`=as.numeric(PB_DI_FI_shannon_phe2$`HIVRNA load`)

PB_DI_FI_shannon_phe_result2=correlation_4index(PB_DI_FI_shannon_phe2)

PB_DI_FI_shannon_phe_result=merge(PB_DI_FI_shannon_phe_result1,PB_DI_FI_shannon_phe_result2,
                                  by=c("x","y"),
                                  all=T)
####correlation between beta diversity and HIV-related phenotypes#### 
for (i in (14:16)){df[,i]=as.factor(df[,i])}

results1=data.frame(factor=NA,R2=NA,p.value=NA)
for(i in c(1:10,14:16,19:21)){
  df1=pheno2[,i,drop=F]
  name=row.names(df1)[which(!is.na(df1[,1]))]
  HIVbray_noNA=HIVbray[row.names(HIVbray)%in%name,colnames(HIVbray)%in%name]
  pp=pheno2[row.names(HIVbray_noNA),]
  PERMANOVA_withinHIV=adonis(HIVbray_noNA~GEN+AGE+log10_read+pp[,i],pp,permutations = 1000)
  PERMANOVA_withinHIV_result=as.data.frame(PERMANOVA_withinHIV$aov.tab)
  results1[i,1]=colnames(df1)[1]
  results1[i,2]=PERMANOVA_withinHIV_result[4,5]
  results1[i,3]=PERMANOVA_withinHIV_result[4,6]
}

results2=data.frame(factor=NA,R2=NA,p.value=NA)
for(i in c(1:10,14:16)){
  df1=pheno2[,i,drop=F]
  name=row.names(df1)[which(!is.na(df1[,1]))]
  name2=row.names(pheno2)[which(!is.na(pheno2$SO)&!is.na(pheno2$`Num-P`))]
  name3=intersect(name,name2)
  HIVbray_noNA=HIVbray[row.names(HIVbray)%in%name3,colnames(HIVbray)%in%name3]
  pp=pheno2[row.names(HIVbray_noNA),]
  PERMANOVA_withinHIV=adonis(HIVbray_noNA~GEN+AGE+log10_read+SO+`Num-P`+pp[,i],pp,permutations = 1000)
  PERMANOVA_withinHIV_result=as.data.frame(PERMANOVA_withinHIV$aov.tab)
  results2[i,1]=colnames(df1)[1]
  results2[i,2]=PERMANOVA_withinHIV_result[6,5]
  results2[i,3]=PERMANOVA_withinHIV_result[6,6]
}
results1$BH.FDR=p.adjust(results1[,'p.value'],method = "BH")
results2$BH.FDR=p.adjust(results2[,'p.value'],method = "BH")
beta_results=merge(results1,results2,
                   by.x = "factor",
                   by.y = "factor",
                   all=T)
beta_results$x=c("beta diversity")
beta_results=beta_results[c(8,1,3,4,2,6,7,5)]
colnames(beta_results)=colnames(PB_DI_FI_shannon_phe_result)
####File: Association between HIV-related phenotypes with alpha and beta diversity, P/B ratio, DI score and FI score####
allresult=rbind(PB_DI_FI_shannon_phe_result,beta_results)
allresult$BH_FDR.x=p.adjust(allresult[,'p.value.x'],method = "BH")
allresult$BH_FDR.y=p.adjust(allresult[,'p.value.y'],method = "BH")
####correlation between individual species/pathways with HIV-related phenotypes#### 
#76species
HIV_noantibio_clr_76_inv=as.data.frame(apply(HIV_noantibio_clr[row.names(HIV_noantibio_clr)%in%DA2_re0.05$microbe,],1,invrank))
for (i in (14:16)){df[,i]=as.numeric(df[,i])}
species_phe=merge(HIV_noantibio_clr_76_inv,df,
                by.x = "row.names",
                by.y = "row.names",
                all=F)
row.names(species_phe)=species_phe$Row.names
species_phe$Row.names=NULL
correlation_sp=function(x){
  results=data.frame(x=NA,y=NA,
                     p.value=NA,BH_FDR=NA,rho=NA)
  k=0
  for(j in c(77:ncol(x))){
    for(i in c(1:76)){
      model=cor.test(x[,i], x[,j],method = "spearman")
      results[(i+k),1]=colnames(x)[i]
      results[(i+k),2]=colnames(x)[j]
      results[(i+k),3]=model$p.value ## p value second row, 4th col
      results[(i+k),5]=model$estimate ## beta coef of DI 
    }
    k=k+80
  }
  results[,'BH_FDR']= p.adjust(results[,'p.value'],method = "BH")
  results<-na.omit(results)
}
DA76sp=correlation_sp(species_phe)

species_phe2=merge(HIV_noantibio_clr_76_inv,df2,
                  by.x = "row.names",
                  by.y = "row.names",
                  all=F)
row.names(species_phe2)=species_phe2$Row.names
species_phe2$Row.names=NULL
DA76sp2=correlation_sp(species_phe2)
DA76sp_result=merge(DA76sp,DA76sp2,
                    by=c("x","y"),
                    all=T)
#163pathways
HIV_pathway_clr_163_inv=as.data.frame(apply(HH2[row.names(HH2)%in%DApathway_result0.05$microbe,],1,invrank))
pathway_phe=merge(HIV_pathway_clr_163_inv,df,
                  by.x = "row.names",
                  by.y = "row.names",
                  all=F)
row.names(pathway_phe)=pathway_phe$Row.names
pathway_phe$Row.names=NULL

correlation_path=function(x){
  results=data.frame(x=NA,y=NA,
                     p.value=NA,BH_FDR=NA,rho=NA)
  k=0
  for(j in c(164:ncol(x))){
    for(i in c(1:163)){
      model=cor.test(x[,i], x[,j],method = "spearman")
      results[(i+k),1]=colnames(x)[i]
      results[(i+k),2]=colnames(x)[j]
      results[(i+k),3]=model$p.value ## p value second row, 4th col
      results[(i+k),5]=model$estimate ## beta coef of DI 
    }
    k=k+200
  }
  results[,'BH_FDR']= p.adjust(results[,'p.value'],method = "BH")
  results<-na.omit(results)
}
DA163path=correlation_path(pathway_phe)

pathway_phe2=merge(HIV_pathway_clr_163_inv,df2,
                   by.x = "row.names",
                   by.y = "row.names",
                   all=F)
row.names(pathway_phe2)=pathway_phe2$Row.names
pathway_phe2$Row.names=NULL
DA163path2=correlation_path(pathway_phe2)
DA163path_result=merge(DA163path,DA163path2,
                    by=c("x","y"),
                    all=T)
####draw####
library(patchwork)

p1=ggplot(data = species_phe2, mapping = aes(y =`CA-RNA/CA-DNA`, x =s__Firmicutes_bacterium_CAG_95))+
  geom_point(aes(y =`CA-RNA/CA-DNA`, x =s__Firmicutes_bacterium_CAG_95))+
  geom_smooth(data =species_phe2, mapping = aes(y =`CA-RNA/CA-DNA`, x = s__Firmicutes_bacterium_CAG_95),method = "lm", formula = y ~ x)
p2=ggplot(data = species_phe2, mapping = aes(x =s__Firmicutes_bacterium_CAG_95, y =`CA-DNA`))+
  geom_point(aes(y =`CA-DNA`, x =s__Firmicutes_bacterium_CAG_95))+
  geom_smooth(data =species_phe2, mapping = aes(x =s__Firmicutes_bacterium_CAG_95, y = `CA-DNA`),method = "lm", formula = y ~ x)
p3=ggplot(data = species_phe2, mapping = aes(y =`CA-DNA`, x =s__Prevotella_sp_CAG_5226 ))+
  geom_point(aes(y =`CA-DNA`, x =s__Prevotella_sp_CAG_5226))+
  geom_smooth(data =species_phe2, mapping = aes(y =`CA-DNA`, x = s__Prevotella_sp_CAG_5226),method = "lm", formula = y ~ x)
p1+p2+p3
####
DA_all_cor=na.omit(allresult[,1:5])
DA_all_cor$estimate=DA_all_cor$rho.x
DA_all_cor$estimate[DA_all_cor$p.value.x>0.05]=0
DA_all_cor$estimate[DA_all_cor$x==c("beta diversity")]=10*DA_all_cor$estimate[DA_all_cor$x==c("beta diversity")]
DA_all_cor$BH_FDR.x[DA_all_cor$BH_FDR.x<0.1]=c("*")
DA_all_cor$BH_FDR.x[!DA_all_cor$BH_FDR.x==c("*")]=NA

DA_all_cor$y=factor(DA_all_cor$y,levels=c("CA-RNA","CA-DNA","CA-RNA/CA-DNA",
                                          "CD4-nadir","CD4 counts","CD4/CD8","CD4-recovery-abs","CD4-recovery-re",
                                          "HIVRNA-zenith","HIVRNA load",
                                          "HIV-duration","Time-diagnosis-cART","cART-duration",
                                          "SO","Num-P","RAI"))
DA_all_cor$x=factor(DA_all_cor$x,levels=c("shannon","beta diversity","log10(P/B_ratio+min(matrix)/2)",
                                          "Dysbiosis_index","Function_index"))

p1=ggplot(DA_all_cor, aes(x=x,y=y, fill=estimate))+ geom_tile(aes(fill = estimate), colour = "transparent")+
  scale_fill_gradient2(low = "#023047",mid="white", high = "#ae2012", 
                       breaks = c(min(DA_all_cor$estimate),0,max(DA_all_cor$estimate)))+
  theme(axis.text.x = element_text(size=12,angle =90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=12))+
  geom_text(data=DA_all_cor,aes(x=x,y=y,label=BH_FDR.x))

DA_all_cor=na.omit(allresult[,c(1,2,6:8)])
DA_all_cor$estimate=DA_all_cor$rho.y
DA_all_cor$estimate[DA_all_cor$p.value.y>0.05]=0
DA_all_cor$BH_FDR.y[DA_all_cor$BH_FDR.y<0.1]=c("*")
DA_all_cor$BH_FDR.y[!DA_all_cor$BH_FDR.y==c("*")]=NA

DA_all_cor$y=factor(DA_all_cor$y,levels=c("CA-RNA","CA-DNA","CA-RNA/CA-DNA",
                                          "CD4-nadir","CD4 counts","CD4/CD8","CD4-recovery-abs","CD4-recovery-re",
                                          "HIVRNA-zenith","HIVRNA load",
                                          "HIV-duration","Time-diagnosis-cART","cART-duration",
                                          "SO","Num-P","RAI"))
DA_all_cor$x=factor(DA_all_cor$x,levels=c("shannon","beta diversity","log10(P/B_ratio+min(matrix)/2)",
                                          "Dysbiosis_index","Function_index"))

p2=ggplot(DA_all_cor, aes(x=x,y=y, fill=estimate))+ geom_tile(aes(fill = estimate), colour = "transparent")+
  scale_fill_gradient2(low = "#023047",mid="white", high = "#ae2012", 
                       breaks = c(min(DA_all_cor$estimate),0,max(DA_all_cor$estimate)))+
  theme(axis.text.x = element_text(size=12,angle =90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=12))+
  geom_text(data=DA_all_cor,aes(x=x,y=y,label=BH_FDR.y))
p1+p2
####draw all index####
p1=ggplot(all_index_3cohort,aes(x=Dysbiosis_index,fill=ID.x,alpha = 1/10))+
  geom_density()
p2=ggplot(all_index_3cohort,aes(x=`log10(P/B_ratio+min(matrix)/2)`,fill=ID.x,alpha = 1/10))+
  geom_density()
p3=ggplot(all_index_3cohort,aes(x=Function_index,fill=ID.x,alpha = 1/10))+
  geom_density()
p1+p2+p3

####P. copri strain associated with HIV-related parameters####
####data:P. copri strain relative abundance and phenotypes corrected for age, sex, reads, SO, and Num-P####
test=merge(cut_ja,df2,
           by.x = "Row.names",
           by.y = "row.names",
           all=F)
test=merge(test,strain_IL10[,c(1,13)],
           by.x = "Row.names",
           by.y = "Row.names",
           all=F)
test$Strain=as.factor(test$Strain)
test1=test[test$Strain==1,]
test2=test[test$Strain==2,]
####lm function to test association between P. copri strain and HIV-related parameters####
linear_regression_strain=function(data){
  results=data.frame(x=NA,y=NA,estimate=NA,SE=NA,lower=NA,upper=NA,p.value=NA)
  for (i in c(12:24)){
    model=lm(data[,i]~data[,25])
    sum1=summary(model)
    results[i,1]=c("P.copri")
    results[i,2]=colnames(data)[i]
    results[i,3]=sum1$coefficients[2,1] ## beta coef of y- varialbe
    results[i,4]=sum1$coefficients[2,2] ## beta coef of y- varialbe
    results[i,5]=confint(model,level=0.95)[2,1]#lower
    results[i,6]=confint(model,level=0.95)[2,2]#upper
    results[i,7]=sum1$coefficients[2,4] ## p value second row, 4th col
  }
  results[,'BH.FDR']=p.adjust(results[,'p.value'],method="BH")
  return(results)
}
#Strain1=control-related
strain1_pheno=linear_regression_strain(test1)
strain2_pheno=linear_regression_strain(test2)
strain_pheno=na.omit(merge(strain1_pheno,strain2_pheno,
                           by.x = "y",
                           by.y = "y",
                           all=F))
####heterogeneity test between 2 strains####
heterogeneity_strain2=function(data){
  results=data.frame(cytokine=NA,heterogeneity.p.value=NA)
  for (i in 1:nrow(data)){
    studlab=c("1","2")
    TE=c(data[i,3],data[i,10])
    seTE=c(data[i,4],data[i,11])
    inDf= data.frame(studlab, TE, seTE)
    heterogeneity_strain_1and2=metagen(TE = TE,
                                       seTE = seTE,
                                       studlab = studlab,
                                       data = inDf,
                                       title = "1 vs.2")
    results[i,1]=data$y[i]
    results[i,2]=heterogeneity_strain_1and2$pval.Q
  }
  return(results)
}
strain_pheno_heterogeneity=heterogeneity_strain2(strain_pheno)
strain_pheno_heterogeneity=merge(strain_pheno,strain_pheno_heterogeneity,
                                 by.x = "y",
                                 by.y = "cytokine",
                                 all=F)




####draw####
p1=ggplot(data = test, mapping = aes(x =s__Prevotella_copri, y =`CD4-recovery-abs`))+
  geom_point(aes(color=Strain))+
  geom_smooth(data =test, mapping = aes(x =s__Prevotella_copri, y = `CD4-recovery-abs`,fill=Strain,color=Strain),method = "lm", formula = y ~ x)+   scale_fill_manual(values=c(`1`= "#3d52a0", `2` = "#d8676e"))+   scale_color_manual(values=c(`1`= "#3d52a0", `2` = "#d8676e"))

p2=ggplot(data = test, mapping = aes(x =s__Prevotella_copri, y =`CD4 counts`))+
  geom_point(aes(color=Strain))+
  geom_smooth(data =test, mapping = aes(x =s__Prevotella_copri, y = `CD4 counts`,fill=Strain,color=Strain),method = "lm", formula = y ~ x)+   scale_fill_manual(values=c(`1`= "#3d52a0", `2` = "#d8676e"))+   scale_color_manual(values=c(`1`= "#3d52a0", `2` = "#d8676e"))

library(patchwork)
p1+p2













