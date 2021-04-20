setwd("~/Desktop/P3 miR-9 Follow up/Fang_Bartel_MCell_2015/")

mir16 <- read.table(file="mir16-variants-energy-cleavage.txt", header = TRUE, sep =  "\t")
mir16_seq <- "GUCAGCAGUGCCUUAGCAGCACGUAAAUAUUGGCGUUAAGAUUCUAAAAUUAUCUCCAGUAUUAACUGUGCUGCUGAAGUAAGGUUGAC"
mir16_seq <- gsub("U", "T",mir16_seq)
mir16_ref <- mir16[which(mir16$variant==mir16_seq),]
#ref_structure <- as.character(mir16_ref$STRUCTURE)
#mir16_ref <- mir16[which(mir16$STRUCTURE==ref_structure),]

mir16$count <- 1
structures_mir16 <- aggregate(mir16$count, by=list(STRUCTURE=mir16$STRUCTURE), FUN=sum)
colnames(structures_mir16) <- c("STRUCTURE","N_structures")

##Loop1 - variant2 
# test_structure <- "((((((......((((((((((((.((((((((.(.(((..........))).).)))))))).)).))))))))))......))))))"
# mir16_test <- mir16[which(mir16$STRUCTURE==test_structure),]
# #write.table(mir16_test, "mir16_test.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# 
# # test_sequences <- "GTCAGCAGXXXCTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAAZZZAGGTTGAC"
# part1 <- "GTCAGCAG"
# # XXX <- "XXX"
# part2 <- "CTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAA"
# # ZZZ <- "ZZZ"
# part3 <- "AGGTTGAC"
# # test_sequences2 <- paste0(part1,XXX,part2,ZZZ,part3)
# # test_sequences
# # test_sequences2==test_sequences
# a <- c("A","T","G","C")
# empty_table <- expand.grid(a,a,a,a,a,a)
# rm(a)
# empty_table$XXX <- paste0(empty_table$Var1,empty_table$Var2,empty_table$Var3)
# empty_table$ZZZ <- paste0(empty_table$Var4,empty_table$Var5,empty_table$Var6)
# 
# empty_table$pos1 <- NA
# empty_table$pos2 <- NA
# empty_table$pos3 <- NA
# i <- 2
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var1[i]=="A") {
#     if (empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }else{
#       empty_table$pos1[i] <- FALSE
#     }}
# 
# 
#   if (empty_table$Var1[i]=="T") {
#     if (empty_table$Var6[i]=="A" | empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
# 
#   if (empty_table$Var1[i]=="G") {
#     if (empty_table$Var6[i]=="C" | empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
# 
#   if (empty_table$Var1[i]=="C") {
#     if (empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var2[i]=="A") {
#     if (empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }else{
#       empty_table$pos2[i] <- FALSE
#     }}
# 
# 
#   if (empty_table$Var2[i]=="T") {
#     if (empty_table$Var5[i]=="A" | empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
# 
#   if (empty_table$Var2[i]=="G") {
#     if (empty_table$Var5[i]=="C" | empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
# 
#   if (empty_table$Var2[i]=="C") {
#     if (empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var3[i]=="A") {
#     if (empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }else{
#       empty_table$pos3[i] <- FALSE
#     }}
# 
# 
#   if (empty_table$Var3[i]=="T") {
#     if (empty_table$Var4[i]=="A" | empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
# 
#   if (empty_table$Var3[i]=="G") {
#     if (empty_table$Var4[i]=="C" | empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
# 
#   if (empty_table$Var3[i]=="C") {
#     if (empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   print(i)
# }
# 
# empty_table$count <- 1
# check1 <- aggregate(empty_table$count, by=list(Var1=empty_table$Var1,Var6=empty_table$Var6,pos1=empty_table$pos1), FUN=sum)
# check2 <- aggregate(empty_table$count, by=list(Var2=empty_table$Var2,Var5=empty_table$Var5,pos2=empty_table$pos2), FUN=sum)
# check3 <- aggregate(empty_table$count, by=list(Var3=empty_table$Var3,Var4=empty_table$Var4,pos3=empty_table$pos3), FUN=sum)
# rm(check1,check2,check3)
# empty_table$variant <- paste0(part1,empty_table$XXX,part2,empty_table$ZZZ,part3)
# 
# empty_table <- empty_table[,c(7:11,13)]
# mir16_test <- merge(mir16_test, empty_table, by="variant")
# zero_pairs <- mir16_test[which(mir16_test$pos1==FALSE & mir16_test$pos2==FALSE & mir16_test$pos3==FALSE),]
# all_pairs <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_1pair_1 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==F),]
# only_1pair_2 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_1pair_3 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_1pair <- rbind(only_1pair_1, only_1pair_2, only_1pair_3)
# rm(only_1pair_1, only_1pair_2, only_1pair_3)
# 
# only_2pair_12 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_2pair_13 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_2pair_23 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_2pair <- rbind(only_2pair_12, only_2pair_13, only_2pair_23)
# rm(only_2pair_12, only_2pair_13, only_2pair_23)
# 
# 
# 
# 
# min1_pair <- mir16_test[!(mir16_test$variant %in% zero_pairs$variant),]
# mean(zero_pairs$ratio)
# mean(only_1pair$ratio)
# mean(only_2pair$ratio)
# mean(all_pairs$ratio)
# 
# wilcox.test(zero_pairs$ratio,min1_pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_1pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_2pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,all_pairs$ratio,alternative = "greater")
# 
# write.table(zero_pairs, "mir16-zero_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_1pair, "mir16-only_1pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_2pair, "mir16-only_2pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(all_pairs, "mir16-all_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# 

##Loop2 - variant 9
# test_structure <- "((((....(((.((((((((((((.((((((((.(.(((..........))).).)))))))).)).)))))))))).)))....))))"
# mir16_test <- mir16[which(mir16$STRUCTURE==test_structure),]
# #write.table(mir16_test, "mir16_test.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# 
# # test_sequences <- "GTCAXXXGTGCCTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAAGTAAZZZTGAC"
# #test_structures <- "((((....(((.((((((((((((.((((((((.(.(((..........))).).)))))))).)).)))))))))).)))....))))"
# 
# part1 <- "GTCA"
# # XXX <- "XXX"
# part2 <- "GTGCCTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAAGTAA"
# # ZZZ <- "ZZZ"
# part3 <- "TGAC"
# # test_sequences2 <- paste0(part1,XXX,part2,ZZZ,part3)
# # test_sequences
# # test_sequences2==test_sequences
# a <- c("A","T","G","C")
# empty_table <- expand.grid(a,a,a,a,a,a)
# rm(a)
# empty_table$XXX <- paste0(empty_table$Var1,empty_table$Var2,empty_table$Var3)
# empty_table$ZZZ <- paste0(empty_table$Var4,empty_table$Var5,empty_table$Var6)
# 
# empty_table$pos1 <- NA
# empty_table$pos2 <- NA
# empty_table$pos3 <- NA
# i <- 2
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var1[i]=="A") {
#     if (empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   
#   
#   if (empty_table$Var1[i]=="T") {
#     if (empty_table$Var6[i]=="A" | empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   
#   if (empty_table$Var1[i]=="G") {
#     if (empty_table$Var6[i]=="C" | empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   
#   if (empty_table$Var1[i]=="C") {
#     if (empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var2[i]=="A") {
#     if (empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   
#   
#   if (empty_table$Var2[i]=="T") {
#     if (empty_table$Var5[i]=="A" | empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   
#   if (empty_table$Var2[i]=="G") {
#     if (empty_table$Var5[i]=="C" | empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   
#   if (empty_table$Var2[i]=="C") {
#     if (empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var3[i]=="A") {
#     if (empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   
#   
#   if (empty_table$Var3[i]=="T") {
#     if (empty_table$Var4[i]=="A" | empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   
#   if (empty_table$Var3[i]=="G") {
#     if (empty_table$Var4[i]=="C" | empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   
#   if (empty_table$Var3[i]=="C") {
#     if (empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   print(i)
# }
# 
# empty_table$count <- 1
# check1 <- aggregate(empty_table$count, by=list(Var1=empty_table$Var1,Var6=empty_table$Var6,pos1=empty_table$pos1), FUN=sum)
# check2 <- aggregate(empty_table$count, by=list(Var2=empty_table$Var2,Var5=empty_table$Var5,pos2=empty_table$pos2), FUN=sum)
# check3 <- aggregate(empty_table$count, by=list(Var3=empty_table$Var3,Var4=empty_table$Var4,pos3=empty_table$pos3), FUN=sum)
# rm(check1,check2,check3)
# empty_table$variant <- paste0(part1,empty_table$XXX,part2,empty_table$ZZZ,part3)
# 
# empty_table <- empty_table[,c(7:11,13)]
# mir16_test <- merge(mir16_test, empty_table, by="variant")
# zero_pairs <- mir16_test[which(mir16_test$pos1==FALSE & mir16_test$pos2==FALSE & mir16_test$pos3==FALSE),]
# all_pairs <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_1pair_1 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==F),]
# only_1pair_2 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_1pair_3 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_1pair <- rbind(only_1pair_1, only_1pair_2, only_1pair_3)
# rm(only_1pair_1, only_1pair_2, only_1pair_3)
# 
# only_2pair_12 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_2pair_13 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_2pair_23 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_2pair <- rbind(only_2pair_12, only_2pair_13, only_2pair_23)
# rm(only_2pair_12, only_2pair_13, only_2pair_23)
# 
# 
# 
# 
# min1_pair <- mir16_test[!(mir16_test$variant %in% zero_pairs$variant),]
# mean(zero_pairs$ratio)
# mean(only_1pair$ratio)
# mean(only_2pair$ratio)
# mean(all_pairs$ratio)
# 
# wilcox.test(zero_pairs$ratio,min1_pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_1pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_2pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,all_pairs$ratio,alternative = "greater")
# 
# write.table(zero_pairs, "mir16-zero_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_1pair, "mir16-only_1pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_2pair, "mir16-only_2pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(all_pairs, "mir16-all_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)

# ##Loop3 - variant 4
# test_structure <- "...(((..(((.((((((((((((.((((((((.(.(((..........))).).)))))))).)).)))))))))).)))..)))..."
# mir16_test <- mir16[which(mir16$STRUCTURE==test_structure),]
# #write.table(mir16_test, "mir16_test.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# 
# # test_sequences <- "XXXAGCAGTGCCTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAAGTAAGGTTZZZ"
# #test_structures <- "...(((..(((.((((((((((((.((((((((.(.(((..........))).).)))))))).)).)))))))))).)))..)))..."
# 
# part1 <- ""
# # XXX <- "XXX"
# part2 <- "AGCAGTGCCTTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGAAGTAAGGTT"
# # ZZZ <- "ZZZ"
# part3 <- ""
# # test_sequences2 <- paste0(part1,XXX,part2,ZZZ,part3)
# # test_sequences
# # test_sequences2==test_sequences
# a <- c("A","T","G","C")
# empty_table <- expand.grid(a,a,a,a,a,a)
# rm(a)
# empty_table$XXX <- paste0(empty_table$Var1,empty_table$Var2,empty_table$Var3)
# empty_table$ZZZ <- paste0(empty_table$Var4,empty_table$Var5,empty_table$Var6)
# 
# empty_table$pos1 <- NA
# empty_table$pos2 <- NA
# empty_table$pos3 <- NA
# i <- 2
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var1[i]=="A") {
#     if (empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   
#   
#   if (empty_table$Var1[i]=="T") {
#     if (empty_table$Var6[i]=="A" | empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   
#   if (empty_table$Var1[i]=="G") {
#     if (empty_table$Var6[i]=="C" | empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   
#   if (empty_table$Var1[i]=="C") {
#     if (empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var2[i]=="A") {
#     if (empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   
#   
#   if (empty_table$Var2[i]=="T") {
#     if (empty_table$Var5[i]=="A" | empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   
#   if (empty_table$Var2[i]=="G") {
#     if (empty_table$Var5[i]=="C" | empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   
#   if (empty_table$Var2[i]=="C") {
#     if (empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var3[i]=="A") {
#     if (empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   
#   
#   if (empty_table$Var3[i]=="T") {
#     if (empty_table$Var4[i]=="A" | empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   
#   if (empty_table$Var3[i]=="G") {
#     if (empty_table$Var4[i]=="C" | empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   
#   if (empty_table$Var3[i]=="C") {
#     if (empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   print(i)
# }
# 
# empty_table$count <- 1
# check1 <- aggregate(empty_table$count, by=list(Var1=empty_table$Var1,Var6=empty_table$Var6,pos1=empty_table$pos1), FUN=sum)
# check2 <- aggregate(empty_table$count, by=list(Var2=empty_table$Var2,Var5=empty_table$Var5,pos2=empty_table$pos2), FUN=sum)
# check3 <- aggregate(empty_table$count, by=list(Var3=empty_table$Var3,Var4=empty_table$Var4,pos3=empty_table$pos3), FUN=sum)
# rm(check1,check2,check3)
# empty_table$variant <- paste0(part1,empty_table$XXX,part2,empty_table$ZZZ,part3)
# 
# empty_table <- empty_table[,c(7:11,13)]
# mir16_test <- merge(mir16_test, empty_table, by="variant")
# zero_pairs <- mir16_test[which(mir16_test$pos1==FALSE & mir16_test$pos2==FALSE & mir16_test$pos3==FALSE),]
# all_pairs <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_1pair_1 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==F),]
# only_1pair_2 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_1pair_3 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_1pair <- rbind(only_1pair_1, only_1pair_2, only_1pair_3)
# rm(only_1pair_1, only_1pair_2, only_1pair_3)
# 
# only_2pair_12 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_2pair_13 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_2pair_23 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_2pair <- rbind(only_2pair_12, only_2pair_13, only_2pair_23)
# rm(only_2pair_12, only_2pair_13, only_2pair_23)
# 
# 
# 
# 
# min1_pair <- mir16_test[!(mir16_test$variant %in% zero_pairs$variant),]
# mean(zero_pairs$ratio)
# mean(only_1pair$ratio)
# mean(only_2pair$ratio)
# mean(all_pairs$ratio)
# 
# wilcox.test(zero_pairs$ratio,min1_pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_1pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_2pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,all_pairs$ratio,alternative = "greater")
# 
# write.table(zero_pairs, "mir16-zero_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_1pair, "mir16-only_1pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_2pair, "mir16-only_2pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(all_pairs, "mir16-all_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# 

# ##Loop variant 3
# test_structure <- "((((((..((...(((((((((((.((((((((.(.(((..........))).).)))))))).)).)))))))))...))..))))))"
# mir16_test <- mir16[which(mir16$STRUCTURE==test_structure),]
# #write.table(mir16_test, "mir16_test.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
#                     
# #test_sequences <- "GTCAGCAGTGXXXTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGZZZTAAGGTTGAC"
# #test_structures <- "((((((..((...(((((((((((.((((((((.(.(((..........))).).)))))))).)).)))))))))...))..))))))"
# 
# part1 <- "GTCAGCAGTG"
# XXX <- "XXX"
# part2 <- "TAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTG"
# ZZZ <- "ZZZ"
# part3 <- "TAAGGTTGAC"
# test_sequences2 <- paste0(part1,XXX,part2,ZZZ,part3)
# # test_sequences
# # test_sequences2==test_sequences
# a <- c("A","T","G","C")
# empty_table <- expand.grid(a,a,a,a,a,a)
# rm(a)
# empty_table$XXX <- paste0(empty_table$Var1,empty_table$Var2,empty_table$Var3)
# empty_table$ZZZ <- paste0(empty_table$Var4,empty_table$Var5,empty_table$Var6)
# 
# empty_table$pos1 <- NA
# empty_table$pos2 <- NA
# empty_table$pos3 <- NA
# i <- 2
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var1[i]=="A") {
#     if (empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }else{
#       empty_table$pos1[i] <- FALSE
#     }}
# 
# 
#   if (empty_table$Var1[i]=="T") {
#     if (empty_table$Var6[i]=="A" | empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
# 
#   if (empty_table$Var1[i]=="G") {
#     if (empty_table$Var6[i]=="C" | empty_table$Var6[i]=="T") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
# 
#   if (empty_table$Var1[i]=="C") {
#     if (empty_table$Var6[i]=="G") {
#       empty_table$pos1[i] <- TRUE
#     }
#     else{
#       empty_table$pos1[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var2[i]=="A") {
#     if (empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }else{
#       empty_table$pos2[i] <- FALSE
#     }}
# 
# 
#   if (empty_table$Var2[i]=="T") {
#     if (empty_table$Var5[i]=="A" | empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
# 
#   if (empty_table$Var2[i]=="G") {
#     if (empty_table$Var5[i]=="C" | empty_table$Var5[i]=="T") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
# 
#   if (empty_table$Var2[i]=="C") {
#     if (empty_table$Var5[i]=="G") {
#       empty_table$pos2[i] <- TRUE
#     }
#     else{
#       empty_table$pos2[i] <- FALSE
#     }}
#   print(i)
# }
# for (i in 1:nrow(empty_table)) {
#   if (empty_table$Var3[i]=="A") {
#     if (empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }else{
#       empty_table$pos3[i] <- FALSE
#     }}
# 
# 
#   if (empty_table$Var3[i]=="T") {
#     if (empty_table$Var4[i]=="A" | empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
# 
#   if (empty_table$Var3[i]=="G") {
#     if (empty_table$Var4[i]=="C" | empty_table$Var4[i]=="T") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
# 
#   if (empty_table$Var3[i]=="C") {
#     if (empty_table$Var4[i]=="G") {
#       empty_table$pos3[i] <- TRUE
#     }
#     else{
#       empty_table$pos3[i] <- FALSE
#     }}
#   print(i)
# }
# 
# empty_table$count <- 1
# check1 <- aggregate(empty_table$count, by=list(Var1=empty_table$Var1,Var6=empty_table$Var6,pos1=empty_table$pos1), FUN=sum)
# check2 <- aggregate(empty_table$count, by=list(Var2=empty_table$Var2,Var5=empty_table$Var5,pos2=empty_table$pos2), FUN=sum)
# check3 <- aggregate(empty_table$count, by=list(Var3=empty_table$Var3,Var4=empty_table$Var4,pos3=empty_table$pos3), FUN=sum)
# rm(check1,check2,check3)
# empty_table$variant <- paste0(part1,empty_table$XXX,part2,empty_table$ZZZ,part3)
# 
# empty_table <- empty_table[,c(7:11,13)]
# mir16_test <- merge(mir16_test, empty_table, by="variant")
# zero_pairs <- mir16_test[which(mir16_test$pos1==FALSE & mir16_test$pos2==FALSE & mir16_test$pos3==FALSE),]
# all_pairs <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_1pair_1 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==F),]
# only_1pair_2 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_1pair_3 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_1pair <- rbind(only_1pair_1, only_1pair_2, only_1pair_3)
# rm(only_1pair_1, only_1pair_2, only_1pair_3)
# 
# only_2pair_12 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==F),]
# only_2pair_13 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==T),]
# only_2pair_23 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==T),]
# only_2pair <- rbind(only_2pair_12, only_2pair_13, only_2pair_23)
# rm(only_2pair_12, only_2pair_13, only_2pair_23)
# 
# 
# 
# 
# min1_pair <- mir16_test[!(mir16_test$variant %in% zero_pairs$variant),]
# mean(zero_pairs$ratio)
# mean(only_1pair$ratio)
# mean(only_2pair$ratio)
# mean(all_pairs$ratio)
# 
# wilcox.test(zero_pairs$ratio,min1_pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_1pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,only_2pair$ratio,alternative = "greater")
# wilcox.test(zero_pairs$ratio,all_pairs$ratio,alternative = "greater")
# 
# write.table(zero_pairs, "mir16-zero_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_1pair, "mir16-only_1pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(only_2pair, "mir16-only_2pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
# write.table(all_pairs, "mir16-all_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)

##Loop variant 5
test_structure <- "((((((.......(((((((((((.((((((((.(.(((..........))).).)))))))).)).))))))))).......))))))"
mir16_test <- mir16[which(mir16$STRUCTURE==test_structure),]
#write.table(mir16_test, "mir16_test.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)

#test_sequences <- "GTCAGCAGTGXXXTAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTGZZZTAAGGTTGAC"
#test_structures <-"((((((.......(((((((((((.((((((((.(.(((..........))).).)))))))).)).))))))))).......))))))"

part1 <- "GTCAGCAGTG"
XXX <- "XXX"
part2 <- "TAGCAGCACGTAAATATTGGCGTTAAGATTCTAAAATTATCTCCAGTATTAACTGTGCTGCTG"
ZZZ <- "ZZZ"
part3 <- "TAAGGTTGAC"
test_sequences2 <- paste0(part1,XXX,part2,ZZZ,part3)
# test_sequences
# test_sequences2==test_sequences
a <- c("A","T","G","C")
empty_table <- expand.grid(a,a,a,a,a,a)
rm(a)
empty_table$XXX <- paste0(empty_table$Var1,empty_table$Var2,empty_table$Var3)
empty_table$ZZZ <- paste0(empty_table$Var4,empty_table$Var5,empty_table$Var6)

empty_table$pos1 <- NA
empty_table$pos2 <- NA
empty_table$pos3 <- NA
i <- 2
for (i in 1:nrow(empty_table)) {
  if (empty_table$Var1[i]=="A") {
    if (empty_table$Var6[i]=="T") {
      empty_table$pos1[i] <- TRUE
    }else{
      empty_table$pos1[i] <- FALSE
    }}
  
  
  if (empty_table$Var1[i]=="T") {
    if (empty_table$Var6[i]=="A" | empty_table$Var6[i]=="G") {
      empty_table$pos1[i] <- TRUE
    }
    else{
      empty_table$pos1[i] <- FALSE
    }}
  
  if (empty_table$Var1[i]=="G") {
    if (empty_table$Var6[i]=="C" | empty_table$Var6[i]=="T") {
      empty_table$pos1[i] <- TRUE
    }
    else{
      empty_table$pos1[i] <- FALSE
    }}
  
  if (empty_table$Var1[i]=="C") {
    if (empty_table$Var6[i]=="G") {
      empty_table$pos1[i] <- TRUE
    }
    else{
      empty_table$pos1[i] <- FALSE
    }}
  print(i)
}
for (i in 1:nrow(empty_table)) {
  if (empty_table$Var2[i]=="A") {
    if (empty_table$Var5[i]=="T") {
      empty_table$pos2[i] <- TRUE
    }else{
      empty_table$pos2[i] <- FALSE
    }}
  
  
  if (empty_table$Var2[i]=="T") {
    if (empty_table$Var5[i]=="A" | empty_table$Var5[i]=="G") {
      empty_table$pos2[i] <- TRUE
    }
    else{
      empty_table$pos2[i] <- FALSE
    }}
  
  if (empty_table$Var2[i]=="G") {
    if (empty_table$Var5[i]=="C" | empty_table$Var5[i]=="T") {
      empty_table$pos2[i] <- TRUE
    }
    else{
      empty_table$pos2[i] <- FALSE
    }}
  
  if (empty_table$Var2[i]=="C") {
    if (empty_table$Var5[i]=="G") {
      empty_table$pos2[i] <- TRUE
    }
    else{
      empty_table$pos2[i] <- FALSE
    }}
  print(i)
}
for (i in 1:nrow(empty_table)) {
  if (empty_table$Var3[i]=="A") {
    if (empty_table$Var4[i]=="T") {
      empty_table$pos3[i] <- TRUE
    }else{
      empty_table$pos3[i] <- FALSE
    }}
  
  
  if (empty_table$Var3[i]=="T") {
    if (empty_table$Var4[i]=="A" | empty_table$Var4[i]=="G") {
      empty_table$pos3[i] <- TRUE
    }
    else{
      empty_table$pos3[i] <- FALSE
    }}
  
  if (empty_table$Var3[i]=="G") {
    if (empty_table$Var4[i]=="C" | empty_table$Var4[i]=="T") {
      empty_table$pos3[i] <- TRUE
    }
    else{
      empty_table$pos3[i] <- FALSE
    }}
  
  if (empty_table$Var3[i]=="C") {
    if (empty_table$Var4[i]=="G") {
      empty_table$pos3[i] <- TRUE
    }
    else{
      empty_table$pos3[i] <- FALSE
    }}
  print(i)
}

empty_table$count <- 1
check1 <- aggregate(empty_table$count, by=list(Var1=empty_table$Var1,Var6=empty_table$Var6,pos1=empty_table$pos1), FUN=sum)
check2 <- aggregate(empty_table$count, by=list(Var2=empty_table$Var2,Var5=empty_table$Var5,pos2=empty_table$pos2), FUN=sum)
check3 <- aggregate(empty_table$count, by=list(Var3=empty_table$Var3,Var4=empty_table$Var4,pos3=empty_table$pos3), FUN=sum)
rm(check1,check2,check3)
empty_table$variant <- paste0(part1,empty_table$XXX,part2,empty_table$ZZZ,part3)

empty_table <- empty_table[,c(7:11,13)]
mir16_test <- merge(mir16_test, empty_table, by="variant")
zero_pairs <- mir16_test[which(mir16_test$pos1==FALSE & mir16_test$pos2==FALSE & mir16_test$pos3==FALSE),]
all_pairs <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==T),]
only_1pair_1 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==F),]
only_1pair_2 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==F),]
only_1pair_3 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==F & mir16_test$pos3==T),]
only_1pair <- rbind(only_1pair_1, only_1pair_2, only_1pair_3)
rm(only_1pair_1, only_1pair_2, only_1pair_3)

only_2pair_12 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==T & mir16_test$pos3==F),]
only_2pair_13 <- mir16_test[which(mir16_test$pos1==T & mir16_test$pos2==F & mir16_test$pos3==T),]
only_2pair_23 <- mir16_test[which(mir16_test$pos1==F & mir16_test$pos2==T & mir16_test$pos3==T),]
only_2pair <- rbind(only_2pair_12, only_2pair_13, only_2pair_23)
rm(only_2pair_12, only_2pair_13, only_2pair_23)




min1_pair <- mir16_test[!(mir16_test$variant %in% zero_pairs$variant),]
mean(zero_pairs$ratio)
mean(only_1pair$ratio)
mean(only_2pair$ratio)
mean(all_pairs$ratio)

wilcox.test(zero_pairs$ratio,min1_pair$ratio,alternative = "greater")
wilcox.test(zero_pairs$ratio,only_1pair$ratio,alternative = "greater")
wilcox.test(zero_pairs$ratio,only_2pair$ratio,alternative = "greater")
wilcox.test(zero_pairs$ratio,all_pairs$ratio,alternative = "greater")

write.table(zero_pairs, "mir16-zero_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
write.table(only_1pair, "mir16-only_1pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
write.table(only_2pair, "mir16-only_2pair-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)
write.table(all_pairs, "mir16-all_pairs-cleavage.txt", sep="\t", append = FALSE, col.names = TRUE, row.names = FALSE)



