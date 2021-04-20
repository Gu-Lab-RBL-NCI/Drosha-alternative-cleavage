setwd("~/Desktop/P3 miR-9 Follow up/Fang_Bartel_MCell_2015/")

mir16 <- read.table(file="mir30-variants-energy-cleavage.txt", header = TRUE, sep =  "\t")
mir16_seq <- "GACAGTGAGCGACTGTAAACATCCTCGACTGGAAGCTGTGAAGCCACAGATGGGCTTTCAGTCGGATGTTTGCAGCTGCCTACTGCC"
mir16_ref <- mir16[which(mir16$variant==mir16_seq),]
ref_structure <- as.character(mir16_ref$STRUCTURE)
mir16 <- mir16[which(mir16$STRUCTURE==ref_structure),]

library(stringdist)
mir16$dist <- stringdist(mir16$variant,mir16_seq)
mir16$posGU <- 99
#"GACAGTGAGCGACTGTAAACATCCTCGACTGGAAGCTGTGAAGCCACAGATGGGCTTTCAGTCppGGATGTTTGCAGCTGCCTACTGCC"
#"..(((((.(((.((((((((((((..((((((((((((((....))).....)))))))))))pp)))))))))))).))).))))).."

basepairs <- c(1,3:6,9,10,13:15,17:24,27:31,33:36)
n <- 36
mir16_GUs_final <- mir16[c(0),]
for (n in basepairs) {
  mir16_test <- mir16
  mir16_test$variant <- paste0(substr(mir16_test$variant,1,63),"pp",substr(mir16_test$variant,64,87))
  mir16_test$site1 <- substr(mir16_test$variant,n,n)
  mir16_test$site2 <- substr(mir16_test$variant,89-n+1,89-n+1)

  
  seqGU <- mir16_test[which(mir16_test$site1=="G" & mir16_test$site2=="T"),]
  seqUG <- mir16_test[which(mir16_test$site1=="T" & mir16_test$site2=="G"),]
  
  seqGU_refIn <- nrow(seqGU[which(seqGU$dist==0),])
  seqUG_refIn <- nrow(seqGU[which(seqUG$dist==0),])
  if (seqGU_refIn==1) {
    seqGU <- seqGU[c(0),]
  }
  if (seqUG_refIn==1) {
    seqUG <- seqUG[c(0),]
  }
  
  mir16_GUs <- c(seqGU$variant, seqUG$variant)
  mir16_GUs <-  gsub("p","",mir16_GUs)
  mir16_GUs <- mir16_GUs[!(mir16_GUs %in% mir16_ref$variant)]
  mir16_GU <- mir16[(mir16$variant %in% mir16_GUs),]
  if (nrow(mir16_GU)>0) {
    mir16_GU$posGU <- n
  }
  mir16_GUs_final <- rbind(mir16_GUs_final,mir16_GU)
  mir16_GUs_final <- unique(mir16_GUs_final)
  
  rm(mir16_test,seqGU,seqUG,mir16_GUs,mir16_GU,seqGU_refIn, seqUG_refIn)
  print(n)
}


mir16_nonGU <- mir16[!(mir16$variant %in% mir16_GUs_final$variant),]
wilcox.test(mir16_nonGU$ratio,mir16_GUs_final$ratio)[3]
write.table(mir16_nonGU, "mir30_nonGU-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
write.table(mir16_GUs_final, "mir30_GUs_final-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
mean(mir16_nonGU$ratio)
mean(mir16_GUs_final$ratio)

mir16_GUs_final$count <- 1
pos_GU <- aggregate(mir16_GUs_final$ratio, by=list(POSITION=mir16_GUs_final$posGU), FUN=mean)
pos_GU_sd <- aggregate(mir16_GUs_final$ratio, by=list(POSITION=mir16_GUs_final$posGU), FUN=sd)
pos_GU_N <- aggregate(mir16_GUs_final$count, by=list(POSITION=mir16_GUs_final$posGU), FUN=sum)
pos_GU <- merge(pos_GU, pos_GU_sd, by="POSITION")
pos_GU <- merge(pos_GU, pos_GU_N, by="POSITION")
write.table(pos_GU, "mir30a_GUs_position-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)


pos1 <- mir16_GUs_final[which(mir16_GUs_final$posGU==1),]
pos2 <- mir16_GUs_final[which(mir16_GUs_final$posGU==2),]

wilcox.test(mir16_nonGU$ratio,pos1$ratio)
wilcox.test(mir16_nonGU$ratio,pos2$ratio)

write.table(zero_pairs, "mir16-zero_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
write.table(only_1pair, "mir16-only_1pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
write.table(only_2pair, "mir16-only_2pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
write.table(all_pairs, "mir16-all_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)




