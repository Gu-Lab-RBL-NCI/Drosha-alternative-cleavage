setwd("~/Desktop/P3 Follow up/Fang_Bartel_MCell_2015/")
setwd("~/Desktop/P3 Follow up/Fang_Bartel_MCell_2015/txt_files/")

getSecondaryStructure<-function(sequence){
  input_data<-c(paste('>', sep=''), as.character(sequence))
  input_file<-'temp-miRNA-seq.fa'
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-'temp-miRNA-secondary-struct.txt'
  #secondary_structure<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  secondary_structure<-system(paste('RNAfold < ', input_file, sep=''), intern=T)[3]
  secondary_structure<-strsplit(secondary_structure, ' ', fixed=T)[[1]][1]
  return(toString(secondary_structure))
}
getEnergy<-function(sequence){
  input_data<-c(paste('>', sep=''), as.character(sequence))
  input_file<- 'temp-miRNA-seq.fa'
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  output_file<-'temp-miRNA-secondary-struct.txt'
  #output_data<-system(paste('RNAfold --infile=', input_file, sep=''), intern=T)[3]
  output_data<-system(paste('RNAfold < ', input_file, sep=''), intern=T)[3]
  energy<-as.numeric(strsplit(output_data, ' ', fixed=T)[[1]][4])
  energy<-gsub(')', '', toString(energy), fixed=T)
  energy<-gsub('(', '', toString(energy), fixed=T)
  if(is.na(as.numeric(energy))){
    energy<-strsplit(output_data, ' ', fixed=T)[[1]][2]
    energy<-gsub(')', '', toString(energy), fixed=T)
    energy<-gsub('(', '', toString(energy), fixed=T)
  }
  if(is.na(as.numeric(energy))){
    energy<-strsplit(output_data, ' ', fixed=T)[[1]][3]
    energy<-gsub(')', '', toString(energy), fixed=T)
    energy<-gsub('(', '', toString(energy), fixed=T)
  }
  return(toString(energy))
}

mir16 <- "GUCAGCAGUGCCUUAGCAGCACGUAAAUAUUGGCGUUAAGAUUCUAAAAUUAUCUCCAGUAUUAACUGUGCUGCUGAAGUAAGGUUGAC"
mir16 <- gsub("U", "T",mir16 )

#Recalculate Drosha processing score
variants_16 <- read.table(file="GSE67937_16_variant_counts.txt", header = TRUE, fill = FALSE, sep =  "\t")
variants_16$input_count <- variants_16$input_count+1
variants_16$selection_count <- variants_16$selection_count+1
variants_16$ratio <- variants_16$selection_count/variants_16$input_count
reference <- variants_16[which(variants_16$variant==mir16),]
variants_16$score <- log2(variants_16$ratio/reference$ratio)
reference$score <- log2(reference$ratio/reference$ratio)
#rm(reference)

#Analyze cleavage fidelity
cl_miR16 <- read.table(file="GSM1658993_16_barcode_cleavage_fragments.txt", header = TRUE, fill = FALSE, sep =  "\t")
cl_miR16$counts <- 1
head(cl_miR16)

summed <- aggregate(cl_miR16$counts, by=list(barcode=cl_miR16$barcode, cleavage_fragment=cl_miR16$cleavage_fragment), FUN=sum)
colnames(summed) <- c("barcode","cleavage_fragment","totalcounts")
summed$cleavage_fragment <- as.character(summed$cleavage_fragment)
head(summed)

variant_count_miR16 <- read.table(file="GSM1658987_16_dictionary.txt", header = FALSE, fill = FALSE, sep =  "\t")
barcount_miR16 <- read.table(file="GSM1658990_16_input_barcode_counts.txt", header = TRUE, fill = FALSE, sep =  "\t")
colnames(variant_count_miR16) <- c("variant","barcode")
colnames(barcount_miR16) <- c("input_count","barcode")
#variant_count_miR16 <- merge(variant_count_miR16, counts, by="barcode")
summed <- merge(summed, variant_count_miR16, by="barcode")
head(summed)
summed$fidelity <- nchar(summed$cleavage_fragment)
summed$variant <- as.character(summed$variant)
summed$variant <- gsub(" ","",summed$variant)
summed$variant <- gsub("1","",summed$variant)
summed$variant <- gsub("2","",summed$variant)
summed$variant <- gsub("3","",summed$variant)
summed$variant <- gsub("4","",summed$variant)
summed$variant <- gsub("5","",summed$variant)
summed$variant <- gsub("6","",summed$variant)
summed$variant <- gsub("7","",summed$variant)
summed$variant <- gsub("8","",summed$variant)
summed$variant <- gsub("9","",summed$variant)
summed$variant <- gsub("0","",summed$variant)

summed_canonical <- summed[which(summed$fidelity==13),]
summed_noncanonical <- summed[which(summed$fidelity!=13),]

total <- aggregate(summed$totalcounts, by=list(variant=summed$variant), FUN=sum)
colnames(total) <- c("variant","totalcounts")

total_canonical <- aggregate(summed_canonical$totalcounts, by=list(variant=summed_canonical$variant), FUN=sum)
total_alternative <- aggregate(summed_noncanonical$totalcounts, by=list(variant=summed_noncanonical$variant), FUN=sum)
colnames(total_canonical) <- c("variant","counts_can")
colnames(total_alternative) <- c("variant","counts_alts")
head(total_canonical)
head(total_alternative)

total <- merge(total, total_canonical, by="variant", all=TRUE)
total <- merge(total, total_alternative, by="variant", all=TRUE)
head(total)

total[is.na(total)] <- 0
total$palternative <- total$counts_alts/total$totalcounts*100
reference2 <- total[which(total$variant==mir16),]
total$ratio <- total$palternative/reference2$palternative
total$GC <- (1-nchar(gsub("G","",gsub("C","",total$variant)))/nchar(total$variant))*100

total_over100 <- total[which(total$totalcounts>=1),]
total_over100 <- total_over100[order(-total_over100$palternative),] 

total_over100$ENERGY <- NA
total_over100$STRUCTURE <- NA


setwd("~/Desktop/P3 Follow up/Fang_Bartel_MCell_2015/txt_files/temp")

for (i in 1:nrow(total_over100)) {
  total_over100$ENERGY[i] <- getEnergy(total_over100$variant[i])
  total_over100$STRUCTURE[i] <- getSecondaryStructure(total_over100$variant[i])
  print (i/nrow(total_over100)*100)
}

total_over100$ENERGY <- as.numeric(total_over100$ENERGY)
write.table(total_over100, "mir16-variants-energy-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)

top10 <- total_over100[1:(nrow(total_over100)/100),]
top25 <- total_over100[1:(nrow(total_over100)/4),]
low25 <- total_over100[(nrow(total_over100)*3/4):nrow(total_over100),]
low10 <- total_over100[(nrow(total_over100)*99/100):nrow(total_over100),]


# do weighted mean and sd to give more weight to variants with more reads!
mean(top25$GC)
sd(top25$GC)
mean(low25$GC)
sd(low25$GC)


total_over100$wt <- total_over100$totalcounts/sum(total_over100$totalcounts)*100

mean(top25$ENERGY)
sd(top25$ENERGY)
nrow(top25)
mean(top10$ENERGY)
sd(top10$ENERGY)
nrow(top10)


mean(low25$ENERGY)
sd(low25$ENERGY)
nrow(low25)
mean(low10$ENERGY)
sd(low10$ENERGY)
nrow(low10)

summary(top25$ENERGY)
summary(low25$ENERGY)



#total_all <- total_over100
total_over100 <- total_all
total_over100 <- total_over100[which(total_over100$totalcounts>=400),]

#install.packages("hexbin")
#install.packages("RColorBrewer")
library(hexbin)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

hexbinplot(ratio~ENERGY,data=total_over100,xlab="Energy kcal/mol",
           ylab="Alternative cleavage ratio",colramp=rf, aspect = 1)

hexbinplot(ratio~log2(totalcounts),data=total_over100,xlab="log2(Total counts)",
           ylab="Alternative cleavage ratio",colramp=rf, aspect = 1)

plot(total_over100$ENERGY, total_over100$ratio,pch = 16, cex = 1.3, col = "blue")
# Multiple Linear Regression Example 
fit <- lm(ratio ~ ENERGY, data=total_over100)
summary(fit) # show results
abline(fit, col="red")

plot(log2(total_over100$totalcounts), total_over100$ratio)


structures_canonical <- aggregate(total_over100$counts_can, by=list(variant=total_over100$STRUCTURE), FUN=sum)
structures_alternative <- aggregate(total_over100$counts_alts, by=list(variant=total_over100$STRUCTURE), FUN=sum)
colnames(structures_canonical) <- c("STRUCTURE","counts_can")
colnames(structures_alternative) <- c("STRUCTURE","counts_alts")
structures <- merge(structures_canonical, structures_alternative, by="STRUCTURE", all = TRUE)
rm(structures_canonical,structures_alternative)
structures$totalcounts <- structures$counts_can+structures$counts_alts
structures$palternative <- structures$counts_alts/structures$totalcounts*100

structures$chck <- regexpr("\\)\\(", gsub("\\.","",structures$STRUCTURE))
structures <- structures[which(structures$chck==-1),]
structures$nt <- nchar(structures$STRUCTURE)

install.packages("dummies")
library(dummies)

# example data
max_nt <- max(structures$nt)
num_str <- nrow(structures)
nt <- matrix(nrow = num_str, ncol = max_nt)
nt <-  as.data.frame(nt)
i <- 1
a <- 1

for (i in 1:nrow(structures)) {
  structure <- structures$STRUCTURE[i]
  for (a in 1:nchar(structure)) {
    bracket <- substr(structure, a, a)
    if(bracket=="."){
      basepair <- 0
    }else {
      basepair <- 1
    }
    nt[i,a] <- as.numeric(basepair)
  }
  print(i/nrow(structures)*100)
}

structures <- cbind(structures,nt)  

sums <- cbind(c(1:max_nt), NA)
sums <- as.data.frame(sums)
a <- 1
i <- 1
for (a in 1:ncol(nt)) {
  colsumm <- 0
  for (i in 1:nrow(nt)) {
    cell <- as.numeric(nt[i,a])
    colsumm <- cell+colsumm
  }
  colsumm <- colsumm/num_str
  sums$V2[a] <- colsumm
  print(a/ncol(nt)*100)
}

hexamers <- expand.grid(n1 = c("(", ")", "."), n2 = c("(", ")", "."), n3 = c("(", ")", "."), n4 = c("(", ")", "."), 
                        n5 = c("(", ")", "."), n6 = c("(", ")", "."))

hexamers$str6 <-paste0(hexamers$n1,hexamers$n2,hexamers$n3,hexamers$n4,hexamers$n5,hexamers$n6)
hexamers$ratio_avg <- NA
hexamers$ratio_sd <- NA
hexamers$n_num <- NA

i <- 1
for (i in 1:nrow(hexamers)) {
  test <- total_over100
  motif <- hexamers$str6[i]
  motif <- gsub("\\(","\\\\(",motif)
  motif <- gsub("\\)","\\\\)",motif)
  motif <- gsub("\\.","\\\\.",motif)
  test$motif_test <- gregexpr(pattern = motif, text=test$STRUCTURE )
  eval <-  test[which(test$motif_test!="-1"),]
  hexamers$ratio_avg[i] <- mean(eval$ratio)
  hexamers$ratio_sd[i] <- sd(eval$ratio)
  hexamers$n_num[i] <- nrow(eval)
  
  print(i/nrow(hexamers)*100)
}





#write.table(total_over100, "total_over100.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
