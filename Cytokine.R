####data####
#all cytokine data overlapping between HIV and 500FG cohort: all
    #invrank-transformed to make a normal distribution:all_inv
#all raw HIV_cytokines (24 deferentially abundant cytokines+ 8 anti-inflammatory cytokines):HIV_cyto
    #correct age+gender+read cout:df, correct age+gender+read cout+SO+Num-P:df2
#shannon index, P/B ration, DI and FI score: all_index_3cohort
#Deferentially abundant species/pathway:HIV_noantibio_clr_76_inv, HIV_pathway_clr_163_inv
#500FG relative abundance of species (clr-transformed):FGspp_spp3_clr

####correlation between cytokines and phenotypes####
pheno_cyto=merge(pheno[,-c(17:18)],HIV_cyto[,1:32],
                 by.x="row.names",
                 by.y="row.names",
                 all=F)
row.names(pheno_cyto)=pheno_cyto$Row.names
pheno_cyto$Row.names=NULL
for (i in names(pheno_cyto)){pheno_cyto[,i]=as.numeric(pheno_cyto[,i])}
library(psych)
cor_13pheno=psych::corr.test(pheno_cyto,method = 'spearman', adjust = 'none')
library(reshape2)
melt_cor_13pheno=cor_13pheno$`r`
melt_p_13pheno=cor_13pheno$`p`
melt_cor_13pheno=melt(melt_cor_13pheno)
melt_p_13pheno=melt(melt_p_13pheno)
melt_pheno=cbind(melt_cor_13pheno,melt_p_13pheno[,3])
colnames(melt_pheno)[3]=c("r")
colnames(melt_pheno)[4]=c("p.value")
melt_pheno$FDR=p.adjust(melt_pheno$p.value,method = "BH")
melt_pheno$final=NA
melt_pheno$final[which(melt_pheno$FDR<0.1)]=melt_pheno$r[which(melt_pheno$FDR<0.1)]
melt_pheno$final[which(is.na(melt_pheno$final))]=0
melt_pheno=melt_pheno[grep("PBMC",melt_pheno$Var2),]
melt_pheno=melt_pheno[-grep("PBMC",melt_pheno$Var1),]

ggplot(melt_pheno, aes(x=Var1,y=Var2, fill=final))+ geom_tile(aes(fill = final), colour = "transparent")+
  scale_fill_gradient2(low = "#023047",mid="white", high = "#ae2012", 
                       breaks = c(min(melt_pheno$final),0,max(melt_pheno$final)))+
  theme(axis.text.x = element_text(angle =90, vjust=0.5, hjust=1))

####Different cytokine productions between PLHIV and HCs from 500FG####
linear_regression_DA_cyto=function(data,re_abundance){
  results=data.frame(cytokine=NA,HIVCohort_p.value=NA,HIVCohort_BH_FDR=NA,HIVCohort_beta=NA,
                     mean_HIV=NA,mean_healthy=NA,"fold change HIV/healthy"=NA)
  for(i in 1:(ncol(data)-5)){
    healthy=re_abundance[grep("HV",row.names(re_abundance)),i]
    hiv=re_abundance[grep("X",row.names(re_abundance)),i]
    data$Group=as.factor(data$Group)
    data$GEN=as.factor(data$GEN)
    model=lm(data[,i]~Group+GEN+AGE,data = data)
    sum1=summary(model)
    results[i,1]=colnames(data)[i]
    results[i,2]=sum1$coefficients[2,4] ## Cohort_p.value
    results[i,4]=sum1$coefficients[2,1] ## Cohort_estimate
    results[i,5]=mean(as.numeric(hiv),na.rm = T)
    results[i,6]=mean(as.numeric(healthy),na.rm = T)
    ratio=results[i,5]/results[i,6]
    results[i,7]=ratio
  }
  results[,3]= p.adjust(results[,'HIVCohort_p.value'],method = "BH")
  results<-results
}
DA_cyto_result=linear_regression_DA_cyto(all_inv,all)

####association between cytokine production and bacterial index####

####residules of cytokine production corrected for gender+age+read counts####
#invrank
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
HIV_cyto_inv=HIV_cyto
for (i in (names(HIV_cyto_inv)[1:32])){HIV_cyto_inv[,i]=as.numeric(HIV_cyto_inv[,i])}
HIV_cyto_inv[,c(1:32)]=lapply(HIV_cyto[,c(1:32)],invrank)
#correct
residule_pheno1=lapply(HIV_cyto_inv[c(1:32)],function(x){resid(lm(x~GEN+AGE+log10_read,HIV_cyto_inv,na.action=na.exclude))})
df<-data.frame(matrix(1:143,ncol=1))
for(i in names(residule_pheno1)){
  row<-data.frame(residule_pheno1[[i]])
  colnames(row)=i
  df<-as.data.frame(cbind(df,row))
}
df=df[,-c(1)]
#residules of cytokine production corrected for gender+age+read counts+SO+Num-P
#correct
HIV_cyto_inv$SO=as.numeric(as.character(HIV_cyto_inv$SO))
residule_pheno2=lapply(HIV_cyto_inv[c(1:32)],function(x){resid(lm(x~GEN+AGE+log10_read+SO+`Num-P`,HIV_cyto_inv,na.action=na.exclude))})
df2<-data.frame(matrix(1:143,ncol=1))
for(i in names(residule_pheno2)){
  row<-data.frame(residule_pheno2[[i]])
  colnames(row)=i
  df2<-as.data.frame(cbind(df2,row))
}
df2=df2[,-c(1)]

####correlation between shannon index, P/B ration, DI and FI score with cytokine#### 
PB_DI_FI_shannon_cyto=merge(all_index_3cohort[,c(1,2,5,6,7)],df,
                           by.x = "Row.names",
                           by.y = "row.names",
                           all=F)
row.names(PB_DI_FI_shannon_cyto)=PB_DI_FI_shannon_cyto$Row.names
PB_DI_FI_shannon_cyto$Row.names=NULL

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
PB_DI_FI_shannon_cyto_result1=correlation_4index(PB_DI_FI_shannon_cyto)

PB_DI_FI_shannon_cyto2=merge(all_index_3cohort[,c(1,2,5,6,7)],df2,
                            by.x = "Row.names",
                            by.y = "row.names",
                            all=F)
row.names(PB_DI_FI_shannon_cyto2)=PB_DI_FI_shannon_cyto2$Row.names
PB_DI_FI_shannon_cyto2$Row.names=NULL

PB_DI_FI_shannon_cyto_result2=correlation_4index(PB_DI_FI_shannon_cyto2)

PB_DI_FI_shannon_cyto_result=merge(PB_DI_FI_shannon_cyto_result1,PB_DI_FI_shannon_cyto_result2,
                                  by=c("x","y"),
                                  all=T)

####correlation between beta diversity and HIV-related phenotypes#### 

results1=data.frame(factor=NA,R2=NA,p.value=NA)
for(i in c(1:32)){
  df1=HIV_cyto[,i,drop=F]
  name=row.names(df1)[which(!is.na(df1[,1]))]
  HIVbray_noNA=HIVbray[row.names(HIVbray)%in%name,colnames(HIVbray)%in%name]
  pp=HIV_cyto[row.names(HIVbray_noNA),]
  PERMANOVA_withinHIV=adonis(HIVbray_noNA~GEN+AGE+log10_read+pp[,i],pp,permutations = 1000)
  PERMANOVA_withinHIV_result=as.data.frame(PERMANOVA_withinHIV$aov.tab)
  results1[i,1]=colnames(df1)[1]
  results1[i,2]=PERMANOVA_withinHIV_result[4,5]
  results1[i,3]=PERMANOVA_withinHIV_result[4,6]
}

results2=data.frame(factor=NA,R2=NA,p.value=NA)
for(i in c(1:32)){
  df1=HIV_cyto[,i,drop=F]
  name=row.names(df1)[which(!is.na(df1[,1]))]
  name2=row.names(HIV_cyto)[which(!is.na(HIV_cyto$SO)&!is.na(HIV_cyto$`Num-P`))]
  name3=intersect(name,name2)
  HIVbray_noNA=HIVbray[row.names(HIVbray)%in%name3,colnames(HIVbray)%in%name3]
  pp=HIV_cyto[row.names(HIVbray_noNA),]
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
colnames(beta_results)=colnames(PB_DI_FI_shannon_cyto_result)
####File: Association between HIV-related phenotypes with alpha and beta diversity, P/B ratio, DI score and FI score####
allresult=rbind(PB_DI_FI_shannon_cyto_result,beta_results)
allresult$BH_FDR.x=p.adjust(allresult[,'p.value.x'],method = "BH")
allresult$BH_FDR.y=p.adjust(allresult[,'p.value.y'],method = "BH")
####correlation between individual species/pathways with HIV-related phenotypes#### 
#76species
HIV_noantibio_clr_76_inv=as.data.frame(apply(HIV_noantibio_clr[row.names(HIV_noantibio_clr)%in%DA2_re0.05$microbe,],1,invrank))
species_cyto=merge(HIV_noantibio_clr_76_inv,df,
                  by.x = "row.names",
                  by.y = "row.names",
                  all=F)
row.names(species_cyto)=species_cyto$Row.names
species_cyto$Row.names=NULL

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
DA76sp=correlation_sp(species_cyto)

species_cyto2=merge(HIV_noantibio_clr_76_inv,df2,
                   by.x = "row.names",
                   by.y = "row.names",
                   all=F)
row.names(species_cyto2)=species_cyto2$Row.names
species_cyto2$Row.names=NULL
DA76sp2=correlation_sp(species_cyto2)
DA76sp_result=merge(DA76sp,DA76sp2,
                    by=c("x","y"),
                    all=T)
#163pathways
HIV_pathway_clr_163_inv=as.data.frame(apply(HH2[row.names(HH2)%in%DApathway_result0.05$microbe,],1,invrank))
pathway_cyto=merge(HIV_pathway_clr_163_inv,df,
                  by.x = "row.names",
                  by.y = "row.names",
                  all=F)
row.names(pathway_cyto)=pathway_cyto$Row.names
pathway_cyto$Row.names=NULL

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
DA163path=correlation_path(pathway_cyto)

pathway_cyto2=merge(HIV_pathway_clr_163_inv,df2,
                   by.x = "row.names",
                   by.y = "row.names",
                   all=F)
row.names(pathway_cyto2)=pathway_cyto2$Row.names
pathway_cyto2$Row.names=NULL
DA163path2=correlation_path(pathway_cyto2)
DA163path_result=merge(DA163path,DA163path2,
                       by=c("x","y"),
                       all=T)
####meta analysis of cytokine-species association between HIV and 500FG cohort####
#500 FG cytokine-species association
FG_cyto=all[all$Group==c("500FG"),c(which(colnames(all)%in%DA_cyto_result$cytokine[DA_cyto_result$HIVCohort_BH_FDR<0.05]),34:38)]
#invrank
FG_cyto_inv=FG_cyto
for (i in (names(FG_cyto_inv)[1:24])){FG_cyto_inv[,i]=as.numeric(FG_cyto_inv[,i])}
FG_cyto_inv[,c(1:24)]=apply(FG_cyto[,c(1:24)],2,invrank)
#correct
residule_pheno3=lapply(FG_cyto_inv[c(1:24)],function(x){resid(lm(x~GEN+AGE+log10_read,FG_cyto_inv,na.action=na.exclude))})
df3<-data.frame(matrix(1:162,ncol=1))
for(i in names(residule_pheno3)){
  row<-data.frame(residule_pheno3[[i]])
  colnames(row)=i
  df3<-as.data.frame(cbind(df3,row))
}
df3=df3[,-c(1)]
df3[,c(1:24)]=apply(df3[,c(1:24)],2,invrank)

fg76=as.data.frame(FGspp_spp3_clr[row.names(FGspp_spp3_clr)%in%colnames(HIV_noantibio_clr_76_inv),colnames(FGspp_spp3_clr)%in%row.names(df3)])
FGspp_spp3_clr_inv=as.data.frame(apply(fg76,1,invrank))

FG_cyto_sp=merge(FGspp_spp3_clr_inv,df3,
                   by.x = "row.names",
                   by.y = "row.names",
                   all=F)
row.names(FG_cyto_sp)=FG_cyto_sp$Row.names
FG_cyto_sp$Row.names=NULL

DA76sp_result_500FG=correlation_sp(FG_cyto_sp)

####meta####

DA76sp_result_500FG$n=162
DA76sp$n=143
HIV_500FG_meta_input=merge(DA76sp,DA76sp_result_500FG,
                      by=c("x","y"),
                      all=F)
meta_res<-my_batch_meta_lm(HIV_500FG_meta_input,c("HIV","500FG"),c(5,9),c(6,10))
meta_res=merge(DA76sp_result,meta_res,
               by=c("x","y"),
               all=T)
#####meta function####
library(tidyverse)
library(reshape2)
library(meta)
library(metafor)
library(circlize)

library(ggsci)
library(viridis)
library(wesanderson)

## Meta-analysis
my_meta_lm <- function(inVec, study_name, beta_col, n_col) {
  require(meta)
  require(metafor)
  
  #inVec<-cbind_vsv_ba_lm_edge[1,]
  #study_name<-c("LLD", "300OB")
  #beta_col<-c(3,10)
  #se_col<-c(4,11)
  
  study_beta <- inVec[beta_col] %>% as.numeric
  study_n <- inVec[n_col] %>% as.numeric
  
  #study_beta<-c(0.7, 0.65)
  #study_se<-c(0.07, 0.08)
  #stydy_n<-c(1000, 300)
  
  
  inDf <- data.frame(study_name, study_beta, study_n)
  
  m.hksj <- metacor(
    cor=study_beta,
    n=study_n,
    data = inDf,
    studlab = study_name
  )
  
  m.hksj.res <- summary(m.hksj)
  
  return(
    list(
      Meta.I2 = m.hksj.res$I2$TE,
      Meta.hetero.Q = m.hksj$Q,
      Meta.hetero.Q.freedom = m.hksj$df.Q,
      Meta.hetero.p = m.hksj$pval.Q
    )
  )
}

## Batch meta-analysis
my_batch_meta_lm <- function(inDf,study_name,beta_col,n_col,row_var_col = 1,col_var_col = 2) {
  #inDf<-assoc
  #study_name<-c("LLD", "300OB")
  #beta_col<-c(4,10)
  #se_col<-c(5,11)
  #row_var_col<-1
  #col_var_col<-2
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_lm,
    study_name = study_name,
    beta_col = beta_col,
    n_col = n_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 4, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.I2', 'Meta.hetero.Q', "Meta.hetero.Q.freedom", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.hetero.BH.FDR = p.adjust(batch_res_edge$Meta.hetero.p, method = "BH")
  )
  
  return(batch_res_edge)
}

####draw the correclation####
p1=ggplot(data = species_cyto, mapping = aes(x =s__Prevotella_copri, y =`Pam3Cys_PBMC_24h_IL.10`))+
  geom_point()+
  geom_smooth(data =species_cyto, mapping = aes(x =s__Prevotella_copri, y = `Pam3Cys_PBMC_24h_IL.10`),method = "lm", formula = y ~ x)

p2=ggplot(data = species_cyto, mapping = aes(x =s__Bacteroides_vulgatus, y =`Pam3Cys_PBMC_24h_IL.1b`))+
  geom_point()+
  geom_smooth(data =species_cyto, mapping = aes(x =s__Bacteroides_vulgatus, y = `Pam3Cys_PBMC_24h_IL.1b`),method = "lm", formula = y ~ x)
p1|p2
ggplot(data = FG_cyto_sp, mapping = aes(x =s__Bacteroides_vulgatus, y =`Pam3Cys_PBMC_24h_IL.1b`))+
  geom_point()+
  geom_smooth(data =FG_cyto_sp, mapping = aes(x =s__Bacteroides_vulgatus, y = `Pam3Cys_PBMC_24h_IL.1b`),method = "lm", formula = y ~ x)


