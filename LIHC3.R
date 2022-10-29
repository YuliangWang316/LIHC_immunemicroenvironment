library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(EDASeq)
library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(org.Hs.eg.db)
library(tidyverse)
library(GSVA)
library(IDConverter)
library("limma")
setwd("E:/TEST3/31.CIBERSORT/")
request_cancer=c("ACC","CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS")
i<-"LIHC"
cancer_type=paste("TCGA",i,sep="-")
print(cancer_type)
clinical<-GDCquery_clinic(project = cancer_type, type = "clinical")
query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query)
count_matrix<-as.data.frame(assay(expdat))
count_gl<-TCGAanalyze_Normalization(count_matrix, geneInfoHT,method =  'geneLength')
remove(count_matrix,expdat,query)
genename<-rownames(count_gl)
e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
e<-e[!duplicated(e$SYMBOL),]
count_gl<-count_gl[e$ENSEMBL,]
rownames(count_gl)<-e$SYMBOL
remove(e,genename)
count_gl<-as.data.frame(count_gl)

newname_gl<-filter_tcga_barcodes(colnames(count_gl),analyte_target = "RNA")
count_gl_new<-count_gl[,newname_gl]
count_gl_new_new<-count_gl_new[rowMeans(count_gl_new)>0,]

v <-voom(count_gl_new_new, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol_LIHC.txt",sep="\t",quote=F,col.names=F)  
source("TMEimmune31.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol_LIHC.txt", perm=100, QN=TRUE)

remove(results,v,out)
remove(newname_gl,count_gl,count_gl_new)
remove(request_cancer,i,cancer_type)

write.table(count_gl_new_new,file="uniq.symbol_LIHC_estimate.txt",sep="\t",quote=F,col.names=TRUE)
library(limma)
library(estimate)
filterCommonGenes(input.f="uniq.symbol_LIHC_estimate.txt", 
                  output.f="commonGenes_LIHC_estimate.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes_LIHC_estimate.gct",
              output.ds="estimateScore_LIHC_estimate.gct", 
              platform="illumina")


scores=read.table("estimateScore_LIHC_estimate.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores)
write.table(out,file="scores_LIHC_estimate.txt",sep="\t",quote=F,col.names=F)

h<-read.table("Msig.txt",sep = "\t")
rownames(h)<-h$V1
h<-as.data.frame(t(h))
h<-h[-1,]
h<-h[-1,]
h<-as.list(h)
gsva<-gsva(expr = as.matrix(count_gl_new_new),gset.idx.list = h,kcdf="Poisson",parallel.sz=20)
write.table(gsva,file="GSVA_LIHC.txt",sep="\t")
remove(h,gsva)
remove(out,scores)
library(tidyverse)
library(DESeq2)

samplesTP <- TCGAquery_SampleTypes(colnames(count_gl_new_new), typesample = c("TP"))
count_gl_new_new_tumor<-count_gl_new_new[,samplesTP]

cibersort<-read.table("CIBERSORT-Results_LIHC.txt",sep = "\t",header = TRUE,row.names = 1)
cibersort_tumor<-cibersort[samplesTP,]

cibersort_tumor1<-data.frame(rownames(cibersort_tumor),cibersort_tumor[,6])
rownames(cibersort_tumor1)<-rownames(cibersort_tumor)
colnames(cibersort_tumor1)<-c("name","CD4_mem_Res")
median1<-median(cibersort_tumor1$CD4_mem_Res)
hi<-rownames(cibersort_tumor1[cibersort_tumor1[,2] > median1,])
lo<-rownames(cibersort_tumor1[cibersort_tumor1[,2] <= median1,])
count_gl_new_new_tumor_new<-count_gl_new_new_tumor[,c(hi,lo)]

condition<-factor(c(rep("hi",length(hi)),rep("lo",length(lo))),levels = c("lo","hi"))
colData<-data.frame(row.names = colnames(count_gl_new_new_tumor_new),condition)

dds <- DESeqDataSetFromMatrix(count_gl_new_new_tumor_new, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="LIHC_CD4_mem_resting_hi_lo_diff.csv")
remove(cibersort_tumor1,colData,dds,res,condition,median1)
diffExp<-read.csv("LIHC_CD4_mem_resting_hi_lo_diff.csv",header = TRUE,row.names = 1)
diffExp_hi<-diffExp[which(diffExp$padj < 0.05 & diffExp$log2FoldChange>1),]
diffExp_lo<-diffExp[which(diffExp$padj < 0.05 & diffExp$log2FoldChange < -1),]
diffExp_new<-count_gl_new_new_tumor_new[c(rownames(diffExp_hi),rownames(diffExp_lo)),]
diffSig<-diffExp[c(rownames(diffExp_hi),rownames(diffExp_lo)),]
geneNum=50     
diffSig=diffSig[order(as.numeric(as.vector(diffSig$log2FoldChange)),decreasing = TRUE),]
diffGeneName=rownames((diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=diffExp_new[hmGene,]
hmExp=log2(hmExp+1)
Type=c(rep("hi",length(hi)),rep("lo",length(lo)))
names(Type)=colnames(diffExp_new)
Type=as.data.frame(Type)
pdf(file="LIHC_CD4_mem_res.heatmap.pdf",height=8,width=10)
library(pheatmap)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()

gene_list=diffExp[,c("log2FoldChange","padj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
# Spgenes<-colored_point[rownames(colored_point) == "Bhlhe41" | rownames(colored_point) == "Arid3a",]
gene_list$threshold<-as.character(gene_list$threshold)
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC > 0)]<-"Up"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC < 0)]<-"Down"
gene_list<-na.omit(gene_list)
Mycolors<-c("Blue","Black","Red")
library("ggplot2")
pdf("vocano_LIHC_CD4_mem_rest.pdf")

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.4, size=1.75)  + xlim(c(-5, 5)) + ylim(c(0, 30)) +xlab("log2 fold change") + ylab("-log10 p-value")+ theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors)
print(g)
dev.off()
remove(Mycolors,lo,hmGene,hi,geneNum,diffLength,diffGeneName,Type,hmExp,gene_list,g,diffSig,diffExp_new,diffExp_lo,diffExp_hi,diffExp,count_gl_new_new_tumor_new)
remove(colored_point)

cibersort_tumor1<-data.frame(rownames(cibersort_tumor),cibersort_tumor[,4])
rownames(cibersort_tumor1)<-rownames(cibersort_tumor)
colnames(cibersort_tumor1)<-c("name","CD8")
median1<-median(cibersort_tumor1$CD8)
hi<-rownames(cibersort_tumor1[cibersort_tumor1[,2] > median1,])
lo<-rownames(cibersort_tumor1[cibersort_tumor1[,2] <= median1,])
count_gl_new_new_tumor_new<-count_gl_new_new_tumor[,c(hi,lo)]

condition<-factor(c(rep("hi",length(hi)),rep("lo",length(lo))),levels = c("lo","hi"))
colData<-data.frame(row.names = colnames(count_gl_new_new_tumor_new),condition)

dds <- DESeqDataSetFromMatrix(count_gl_new_new_tumor_new, colData, design= ~ condition)
dds <- DESeq(dds)

res= results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="LIHC_CD8_hi_lo_diff.csv")
remove(cibersort_tumor1,colData,dds,res,condition,median1)
diffExp<-read.csv("LIHC_CD8_hi_lo_diff.csv",header = TRUE,row.names = 1)
diffExp_hi<-diffExp[which(diffExp$padj < 0.05 & diffExp$log2FoldChange>1),]
diffExp_lo<-diffExp[which(diffExp$padj < 0.05 & diffExp$log2FoldChange < -1),]
diffExp_new<-count_gl_new_new_tumor_new[c(rownames(diffExp_hi),rownames(diffExp_lo)),]
diffSig<-diffExp[c(rownames(diffExp_hi),rownames(diffExp_lo)),]
geneNum=50     
diffSig=diffSig[order(as.numeric(as.vector(diffSig$log2FoldChange)),decreasing = TRUE),]
diffGeneName=rownames((diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=diffExp_new[hmGene,]
hmExp=log2(hmExp+1)
Type=c(rep("hi",length(hi)),rep("lo",length(lo)))
names(Type)=colnames(diffExp_new)
Type=as.data.frame(Type)
pdf(file="LIHC_CD8.heatmap.pdf",height=8,width=10)
library(pheatmap)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()

gene_list=diffExp[,c("log2FoldChange","padj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
# Spgenes<-colored_point[rownames(colored_point) == "Bhlhe41" | rownames(colored_point) == "Arid3a",]
gene_list$threshold<-as.character(gene_list$threshold)
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC > 0)]<-"Up"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC < 0)]<-"Down"
gene_list<-na.omit(gene_list)
Mycolors<-c("Blue","Black","Red")
library("ggplot2")
pdf("vocano_LIHC_CD8.pdf")

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.4, size=1.75)  + xlim(c(-4, 4)) + ylim(c(0, 40)) +xlab("log2 fold change") + ylab("-log10 p-value")+ theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors)
print(g)
dev.off()
remove(Mycolors,lo,hmGene,hi,geneNum,diffLength,diffGeneName,Type,hmExp,gene_list,g,diffSig,diffExp_new,diffExp_lo,diffExp_hi,diffExp,count_gl_new_new_tumor_new)
remove(colored_point)
Estimate_score<-read.table("scores_LIHC_estimate.txt",sep = "\t",header = TRUE,row.names = 1)
Estimate_score_tumor<-Estimate_score[samplesTP,]
Estimate_score_tumor$sample <- sapply(strsplit(rownames(Estimate_score_tumor),'-'),function(x) paste0(x[1:3],collapse="-"))

clinical$Stromal_Score <- Estimate_score_tumor[match(clinical$submitter_id,Estimate_score_tumor$sample),][,1]
clinical$Immune_Score <- Estimate_score_tumor[match(clinical$submitter_id,Estimate_score_tumor$sample),][,2]
df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,Stromal_Score,age_at_index,gender,ajcc_pathologic_stage))
df <- df[!is.na(df$Stromal_Score),]
df<-df[which(df$vital_status!="NA"),]
for (j in 1:length(rownames(df))) {
  if(is.na(df$days_to_death[j])){
    df$Time[j] <- df$days_to_last_follow_up[j]
  }else if(is.na(df$days_to_last_follow_up[j]) ){
    df$Time[j] <- df$days_to_death[j]
  }
  else if(df$days_to_death[j] >=df$days_to_last_follow_up[j]){
    df$Time[j] <-df$days_to_death[j]
  }
}
df<-df[which(df$Time != 0),]
for (j in 1:length(rownames(df))) {
  if(df$vital_status[j] == "Alive"){
    df$events[j]<-0
  }else if(df$vital_status[j] == "Dead"){
    df$events[j]<-1
  }
}
res.cut<-surv_cutpoint(df,time = "Time",event = "events",variables = "Stromal_Score" )
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(Time,events)~ Stromal_Score,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
res.cox<-coxph(Surv(Time,events) ~ Stromal_Score,data=res.cat)
summary(res.cox)
test.ph<-cox.zph(res.cox)
ggcoxzph(test.ph)
newrescat<-cbind(res.cat,df$age_at_index,df$gender,df$ajcc_pathologic_stage)
# res.cox_new<-coxph(Surv(Time,events) ~ Treg_1C_gc + df$age_at_index + df$gender + df$ajcc_pathologic_stage,data=newrescat)
colnames(newrescat)[4:6]<-c("Age","gender","stage")
res.cox_new<-coxph(Surv(Time,events) ~ Stromal_Score + Age ,data=newrescat)
summary(res.cox_new)
test.ph_new<-cox.zph(res.cox_new)
ggcoxzph(test.ph_new)

ggadjustedcurves(res.cox_new,data = newrescat,variable = "Stromal_Score",risk.table = TRUE,pval = TRUE)

df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,Immune_Score,age_at_index,gender,ajcc_pathologic_stage))
df <- df[!is.na(df$Immune_Score),]
df<-df[which(df$vital_status!="NA"),]
for (j in 1:length(rownames(df))) {
  if(is.na(df$days_to_death[j])){
    df$Time[j] <- df$days_to_last_follow_up[j]
  }else if(is.na(df$days_to_last_follow_up[j]) ){
    df$Time[j] <- df$days_to_death[j]
  }
  else if(df$days_to_death[j] >=df$days_to_last_follow_up[j]){
    df$Time[j] <-df$days_to_death[j]
  }
}
df<-df[which(df$Time != 0),]
for (j in 1:length(rownames(df))) {
  if(df$vital_status[j] == "Alive"){
    df$events[j]<-0
  }else if(df$vital_status[j] == "Dead"){
    df$events[j]<-1
  }
}
res.cut<-surv_cutpoint(df,time = "Time",event = "events",variables = "Immune_Score" )
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(Time,events)~ Immune_Score,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
res.cox<-coxph(Surv(Time,events) ~ Immune_Score,data=res.cat)
summary(res.cox)
test.ph<-cox.zph(res.cox)
ggcoxzph(test.ph)
newrescat<-cbind(res.cat,df$age_at_index,df$gender,df$ajcc_pathologic_stage)
colnames(newrescat)[4:6]<-c("Age","gender","stage")
res.cox_new<-coxph(Surv(Time,events) ~ Immune_Score + Age ,data=newrescat)
summary(res.cox_new)
test.ph_new<-cox.zph(res.cox_new)
ggcoxzph(test.ph_new)

ggadjustedcurves(res.cox_new,data = newrescat,variable = "Immune_Score",risk.table = TRUE,pval = TRUE)

cibersort_tumor$sample <- sapply(strsplit(rownames(cibersort_tumor),'-'),function(x) paste0(x[1:3],collapse="-"))
K <- cibersort_tumor[match(clinical$submitter_id,cibersort_tumor$sample),][,c(1:22)]
clinical<-cbind(clinical,K)
remove(K)
remove(Estimate_score,Estimate_score_tumor)

df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,Neutrophils,age_at_index,gender,ajcc_pathologic_stage))
df <- df[!is.na(df$Neutrophils),]
df<-df[which(df$vital_status!="NA"),]
for (j in 1:length(rownames(df))) {
  if(is.na(df$days_to_death[j])){
    df$Time[j] <- df$days_to_last_follow_up[j]
  }else if(is.na(df$days_to_last_follow_up[j]) ){
    df$Time[j] <- df$days_to_death[j]
  }
  else if(df$days_to_death[j] >=df$days_to_last_follow_up[j]){
    df$Time[j] <-df$days_to_death[j]
  }
}
df<-df[which(df$Time != 0),]
for (j in 1:length(rownames(df))) {
  if(df$vital_status[j] == "Alive"){
    df$events[j]<-0
  }else if(df$vital_status[j] == "Dead"){
    df$events[j]<-1
  }
}
res.cut<-surv_cutpoint(df,time = "Time",event = "events",variables = "Neutrophils" )
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(Time,events)~ Neutrophils,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
res.cox<-coxph(Surv(Time,events) ~ Neutrophils,data=res.cat)
summary(res.cox)
test.ph<-cox.zph(res.cox)
ggcoxzph(test.ph)
newrescat<-cbind(res.cat,df$age_at_index,df$gender,df$ajcc_pathologic_stage)
colnames(newrescat)[4:6]<-c("Age","gender","stage")
res.cox_new<-coxph(Surv(Time,events) ~ Neutrophils + Age ,data=newrescat)
summary(res.cox_new)
test.ph_new<-cox.zph(res.cox_new)
ggcoxzph(test.ph_new)

ggadjustedcurves(res.cox_new,data = newrescat,variable = "Neutrophils",risk.table = TRUE,pval = TRUE)
remove(test.ph,test.ph_new,res.cat,res.cut,res.cox,res.cox_new,newrescat,fit,df)
remove(j)
write.table(clinical,"LIHC_clinical_with_estimate_cibersort.txt",sep = "\t")

library(psych)
cibersort<-read.table("CIBERSORT-Results_LIHC.txt",sep = "\t",header = TRUE,row.names = 1)
cibersort_tumor<-cibersort[samplesTP,]
GSVA_origin<-read.table("GSVA_LIHC.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
GSVA_tumor<-GSVA_origin[,rownames(cibersort_tumor)]

result<-corr.test(as.data.frame(t(count_gl_new_new_tumor)),cibersort_tumor[,1:22],use = "pairwise",method = "spearman",adjust = "fdr",alpha = .05,ci = TRUE,minlength = 10)
write.table(result$r,"LIHC_count_cibersort_corr.txt",sep = "\t")
write.table(result$p.adj,"LIHC_count_cibersort_padj.txt",sep = "\t")
