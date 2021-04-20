library(biomaRt)
MartMir <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#View(listAttributes(MartMir))
ArrayAtt <- getBM(attributes=c("external_gene_name","ensembl_gene_id","ensembl_gene_id_version","ensembl_transcript_id","external_transcript_name","ensembl_transcript_id_version","transcript_biotype","arrayexpress"),mart=MartMir)

setwd("~/Desktop/P3 miR-9 Follow up/TCGA_BRCA")
plotted <- read.table(file="BRCA_miRNA_normal_tumor_SEED_num_v2.txt", header = TRUE, sep =  "\t")
high_fidelity <- plotted[c(1:30),]
low_fidelity <- plotted[c(76:100),]
high_fidelity$group <- "high_fidelity"
low_fidelity$group <- "low_fidelity"
selected_miRNA <- rbind(high_fidelity,low_fidelity)
selected_miRNA <- selected_miRNA[,c(1,9)]
rm(plotted,high_fidelity,low_fidelity)

setwd("~/Desktop/R_scripts/TUMR_targets/Acong-mir21/BRCA-RNAseq-expr/")
BRCA <- read.table(file="BRCA_Normalized_FPKM-UQ.tsv", header = TRUE, sep =  "\t")
metadata <- read.table(file="1516397366532-manifest.csv", header = TRUE, sep =  ",")
metadata <- metadata[,c(2,9,11)]
metadata$name <- gsub(".gz","",metadata$name)


####START HERE###
#################
name_gene <- "SRSF3"
rm(array_test,genes_table,ensg_gene)
array_test <- ArrayAtt
array_test <- array_test[which(array_test$external_gene_name==name_gene),c(1,2)]
array_test <- unique(array_test)
genes_table <- as.data.frame(colnames(BRCA))
genes_table$ensembl_gene_id <- substr(start = 1,stop = 15, genes_table$`colnames(BRCA)`)
genes_table <- merge(genes_table, array_test, by="ensembl_gene_id")
ensg_gene <- as.character(genes_table$`colnames(BRCA)`)


mini <- BRCA[,c("name",ensg_gene)]
mini <- na.omit(mini)
mini <- merge(metadata,mini,by="name")
mini$FPKM <- mini[,4]

normal <- mini[which(mini$sample_type=="Solid Tissue Normal"),]
mean(normal$FPKM)
sd(normal$FPKM)

primary <- mini[which(mini$sample_type=="Primary Tumor"),]
mean(primary$FPKM)
sd(primary$FPKM)
  
setwd("~/Desktop/P3 miR-9 Follow up/TCGA_BRCA/tumor/sequence_tsv/")
# fileNames <- Sys.glob("*.tsv")
# file <- "tumor_BRCA_2_5.isomir.sequence_info.tsv"
# 
# newdata <- read.table(file=file, header = TRUE, sep = "\t")
# newdata$SEED <- substr(newdata$SEQUENCE,2,7)
# MIRNA_seeds <- aggregate(newdata$READS, by=list(MIRNA=newdata$MIRNA, SEED=newdata$SEED), FUN=sum)
# MIRNA_seeds$templated <- NA
# i <- 2
# for (i in 1:nrow(MIRNA_seeds)) {
#   MIRNA <- as.character(MIRNA_seeds$MIRNA[i])
#   SEED <- as.character(MIRNA_seeds$SEED[i])
#   miRBase_test <- miRBase[which(miRBase$MIRNA==MIRNA),]
#   miRBase_test$templated_seed <- regexpr(SEED, miRBase_test$PRI.SEQUENCE)
#   miRBase_test <- miRBase_test[which(miRBase_test$templated_seed!=-1),]
#   MIRNA_seeds$templated[i] <- nrow(miRBase_test)
#   print(i/nrow(MIRNA_seeds)*100)
#   rm(MIRNA,SEED,miRBase_test)
# }
# #MIRNA_seeds <- MIRNA_seeds[which(MIRNA_seeds$templated>0),] #this step filtered 10% of untemplated seeds (in this example)
# 
# MIRNA_counts <- aggregate(MIRNA_seeds$x, by=list(MIRNA=MIRNA_seeds$MIRNA), FUN=sum)
# MIRNA_seeds <- merge(MIRNA_seeds, MIRNA_counts,by="MIRNA")
# MIRNA_seeds$p <- MIRNA_seeds$x.x/MIRNA_seeds$x.y
# MIRNA_seeds$p2 <- MIRNA_seeds$p*MIRNA_seeds$p
# MIRNA_sump2 <- aggregate(MIRNA_seeds$p2, by=list(MIRNA=MIRNA_seeds$MIRNA), FUN=sum)
# MIRNA_sump2$effSEED <- 1/MIRNA_sump2$x
# MIRNA_sump2 <- MIRNA_sump2[,c(1,3)]
# 
# MIRNA_counts$CPM_miRNA <- MIRNA_counts$x/sum(MIRNA_counts$x)*1000000
# MIRNA_counts <- MIRNA_counts[,c(1,3)]
# MIRNA_sump2 <- merge(MIRNA_sump2, MIRNA_counts,by="MIRNA")
# 
# sample_name <- as.character(unique(newdata$SAMPLE))
# sample_name <- gsub("_mirna_gdc_realn.pe_0.fastq","",sample_name)
# MIRNA_sump2$SAMPLE <- "fake"
# MIRNA_seeds_final <- MIRNA_sump2
# rm(file,sample_name, MIRNA_counts, MIRNA_seeds,MIRNA_sump2,newdata,i)
# miRNA_list <- as.character(unique(MIRNA_seeds_final$MIRNA))
# 
# n <- 1
# for (file in fileNames) {
#   newdata <- read.table(file=file, header = TRUE, sep = "\t")
#   newdata <- newdata[(newdata$MIRNA %in% miRNA_list),]
#   newdata$SEED <- substr(newdata$SEQUENCE,2,7)
#   MIRNA_seeds <- aggregate(newdata$READS, by=list(MIRNA=newdata$MIRNA, SEED=newdata$SEED), FUN=sum)
#   MIRNA_seeds$templated <- NA
#   for (i in 1:nrow(MIRNA_seeds)) {
#     MIRNA <- as.character(MIRNA_seeds$MIRNA[i])
#     SEED <- as.character(MIRNA_seeds$SEED[i])
#     miRBase_test <- miRBase[which(miRBase$MIRNA==MIRNA),]
#     miRBase_test$templated_seed <- regexpr(SEED, miRBase_test$PRI.SEQUENCE)
#     miRBase_test <- miRBase_test[which(miRBase_test$templated_seed!=-1),]
#     MIRNA_seeds$templated[i] <- nrow(miRBase_test)
#     rm(MIRNA,SEED,miRBase_test)
#   }
#   MIRNA_seeds <- MIRNA_seeds[which(MIRNA_seeds$templated>0),]
#   MIRNA_counts <- aggregate(MIRNA_seeds$x, by=list(MIRNA=MIRNA_seeds$MIRNA), FUN=sum)
#   
#   MIRNA_seeds <- merge(MIRNA_seeds, MIRNA_counts,by="MIRNA")
#   MIRNA_seeds$p <- MIRNA_seeds$x.x/MIRNA_seeds$x.y
#   MIRNA_seeds$p2 <- MIRNA_seeds$p*MIRNA_seeds$p
#   MIRNA_sump2 <- aggregate(MIRNA_seeds$p2, by=list(MIRNA=MIRNA_seeds$MIRNA), FUN=sum)
#   MIRNA_sump2$effSEED <- 1/MIRNA_sump2$x
#   MIRNA_sump2 <- MIRNA_sump2[,c(1,3)]
#   
#   MIRNA_counts$CPM_miRNA <- MIRNA_counts$x/sum(MIRNA_counts$x)*1000000
#   MIRNA_counts <- MIRNA_counts[,c(1,3)]
#   MIRNA_sump2 <- merge(MIRNA_sump2, MIRNA_counts,by="MIRNA")
#   
#   sample_name <- as.character(unique(newdata$SAMPLE))
#   sample_name <- gsub("_mirna_gdc_realn.pe_0.fastq","",sample_name)
#   MIRNA_sump2$SAMPLE <- sample_name
#   MIRNA_sump2$count <- 1
#   MIRNA_seeds_final$count <- 1
#   MIRNA_seeds_final <- MIRNA_seeds_final[which(MIRNA_seeds_final$SAMPLE!="fake"),]
#   MIRNA_seeds_final <- rbind(MIRNA_sump2,MIRNA_seeds_final)
#   MIRNA_rep <- aggregate(MIRNA_seeds_final$count, by=list(MIRNA=MIRNA_seeds_final$MIRNA), FUN=sum)
#   max_samples <- max(MIRNA_rep$x)
#   MIRNA_rep <- MIRNA_rep[which(MIRNA_rep$x==max_samples),]
#   MIRNA_seeds_final <- MIRNA_seeds_final[(MIRNA_seeds_final$MIRNA %in% MIRNA_rep$MIRNA),]
#   miRNA_list <- as.character(unique(MIRNA_seeds_final$MIRNA))
#   
#   rm(file,sample_name, MIRNA_counts, MIRNA_seeds,MIRNA_sump2,newdata,i)
#   print(n)
#   n <- n+1
# }
setwd("~/Desktop/P3 miR-9 Follow up/TCGA_miRNA_mutants/")
# write.table(MIRNA_seeds_final, "BRCA_seeds_by_sample.txt", sep="\t", append = FALSE, col.names = T, row.names = FALSE)
MIRNA_seeds_final <- read.table(file="BRCA_seeds_by_sample.txt", header = TRUE, sep =  "\t")
MIRNA_seeds_final$sample_id <- substr(MIRNA_seeds_final$SAMPLE, 1, 16)
MIRNA_seeds_final$MIRNA <- gsub("hsa-","",MIRNA_seeds_final$MIRNA)
MIRNA_seeds_final$MIRNA <- gsub("p-1-2-3","p",MIRNA_seeds_final$MIRNA)
MIRNA_seeds_final$MIRNA <- gsub("p-1-2","p",MIRNA_seeds_final$MIRNA)

MIRNA_seeds_final <- MIRNA_seeds_final[MIRNA_seeds_final$sample_id %in% mini$sample_id,]
mini <- mini[mini$sample_id %in% MIRNA_seeds_final$sample_id,]

mini <- mini[order(mini$FPKM),]
n <- 100
m <- nrow(mini)-n+1
Low_Gene <- mini[c(1:n),]
mean(Low_Gene$FPKM)
sd(Low_Gene$FPKM)
High_Gene <- mini[c(m:nrow(mini)),]
mean(High_Gene$FPKM)
sd(High_Gene$FPKM)
Low_Gene <- merge(Low_Gene, MIRNA_seeds_final, by="sample_id")
High_Gene <- merge(High_Gene, MIRNA_seeds_final, by="sample_id")

summary_gene <- aggregate(High_Gene$CPM_miRNA, by=list(MIRNA=High_Gene$MIRNA), FUN=mean)
summary_gene$eff_Low <- NA
summary_gene$eff_High <- NA
summary_gene$effSEED_p <- NA
summary_gene$CPM_Low <- NA
summary_gene$CPM_High <- NA
summary_gene$CPM_p <- NA
i <- 1
for (i in 1:nrow(summary_gene)) {
  miRNA_test <- as.character(summary_gene$MIRNA[i])
  test_low <- Low_Gene[which(Low_Gene$MIRNA==miRNA_test),]
  test_high <- High_Gene[which(High_Gene$MIRNA==miRNA_test),]
  summary_gene$CPM_Low[i] <- mean(test_low$CPM_miRNA)
  summary_gene$CPM_High[i] <- mean(test_high$CPM_miRNA)
  expr_t_p <- t.test(test_low$CPM_miRNA, test_high$CPM_miRNA)[3]
  effSEED_t_p <- t.test(test_low$effSEED, test_high$effSEED)[3]
  summary_gene$eff_Low[i] <- mean(test_low$effSEED)
  summary_gene$eff_High[i] <- mean(test_high$effSEED)
  
  summary_gene$effSEED_p[i] <- as.numeric(effSEED_t_p)
  summary_gene$CPM_p[i] <- as.numeric(expr_t_p)
  print(i)
}
summary_gene$FC_expr <- summary_gene$CPM_High/summary_gene$CPM_Low
summary_gene$log10_expr <- -log10(summary_gene$CPM_p)

summary_gene$FC_NSEED <- summary_gene$eff_High/summary_gene$eff_Low
summary_gene$log10_NSEED <- -log10(summary_gene$effSEED_p)
plot(summary_gene$FC_expr, summary_gene$log10_expr)
plot(summary_gene$FC_NSEED, summary_gene$log10_NSEED, xlim=c(0.8,1.3), ylim=c(0,10))

summary_gene <- merge(summary_gene, selected_miRNA, by="MIRNA", all = T)
high_fidelity <- summary_gene[which(summary_gene$group=="high_fidelity"),]
low_fidelity <- summary_gene[which(summary_gene$group=="low_fidelity"),]
wilcox.test(high_fidelity$FC_expr, low_fidelity$FC_expr, alternative = "less")
wilcox.test(high_fidelity$FC_NSEED, low_fidelity$FC_NSEED, alternative = "greater")
summary_gene <- summary_gene[order(summary_gene$group),]

setwd("~/Desktop/P3 miR-9 Follow up/TCGA_miRNA_mutants/new_factors/")
rm(name_file)
name_file <- paste0(name_gene,"_",ensg_gene,"_volcano.txt")
print(name_file)
write.table(summary_gene, name_file, sep="\t", append = FALSE, col.names = T, row.names = FALSE)

### END






















hist(arm3p$CNNC)
#CNNC <- arm3p[which(arm3p$CNNC>-1),]
CNNC <- arm3p
CNNC <- CNNC[,c(3,26)]
write.table(CNNC, "high_conf_miRNA_flanking_3.txt", sep="\t", append = FALSE, col.names = T, row.names = FALSE)
primiRNA <- miRBase[,c(1,3)]
CNNC <- merge(primiRNA,CNNC,by="PRIMIRNA")



motif_expression <- summary_gene[which(summary_gene$CPM_p<0.05 & summary_gene$FC_expr>1),]
#CPM_SRSF3_expression <- CPM_SRSF3_high[which(CPM_SRSF3_high$CPM_p<0.05),]
#CPM_SRSF3_effSEED <- CPM_SRSF3_high[which(CPM_SRSF3_high$effSEED_p<0.05  & CPM_SRSF3_high$FC_NSEED>1),]
CPM_SRSF3_effSEED <- CPM_SRSF3_high[which(CPM_SRSF3_high$effSEED_p<0.05),]


motif_SRSF3_expression <- merge(CPM_SRSF3_expression,CNNC, by="MIRNA")
motif_SRSF3_expression <- unique(motif_SRSF3_expression[,c(13,14)])
write.table(motif_SRSF3_expression, "motif_SRSF3_expression.txt", sep="\t", append = FALSE, col.names = T, row.names = FALSE)

motif_SRSF3_SEED <- merge(CPM_SRSF3_effSEED,CNNC, by="MIRNA")
motif_SRSF3_SEED <- unique(motif_SRSF3_SEED[,c(13,14)])
write.table(motif_SRSF3_SEED, "motif_SRSF3_NSEED.txt", sep="\t", append = FALSE, col.names = T, row.names = FALSE)




SEED_num_Low_SRSF3 <- aggregate(Low_SRSF3$effSEED, by=list(MIRNA=Low_SRSF3$MIRNA), FUN=mean)
colnames(SEED_num_Low_SRSF3) <- c("MIRNA","effSEED_Low_SRSF3")
SEED_num_High_SRSF3 <- aggregate(High_SRSF3$effSEED, by=list(MIRNA=High_SRSF3$MIRNA), FUN=mean)
colnames(SEED_num_High_SRSF3) <- c("MIRNA","effSEED_High_SRSF3")
SRSF3 <- merge(SEED_num_High_SRSF3, SEED_num_Low_SRSF3, by="MIRNA")
wilcox.test(SRSF3$effSEED_High_SRSF3, SRSF3$effSEED_Low_SRSF3)
SRSF3$delta <- SRSF3$effSEED_Low_SRSF3-SRSF3$effSEED_High_SRSF3

SRSF3_CNNC <- merge(SRSF3, primiRNA, by="MIRNA")
SRSF3_CNNC <- SRSF3_CNNC[SRSF3_CNNC$PRIMIRNA %in% CNNC$PRIMIRNA,]
SRSF3_CNNC <- unique(SRSF3_CNNC[,c(1:4)])
wilcox.test(SRSF3_CNNC$effSEED_High_SRSF3, SRSF3_CNNC$effSEED_Low_SRSF3)
