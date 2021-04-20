setwd("~/Desktop/P3 miR-9 Follow up/mouse_brain/")

cortex1 <- read.table(file="SRR306538_neocortex1.fastq.collapsed.tsv", header = FALSE, fill = FALSE, sep =  "")
cortex1$sample <- "cortex1"
cortex2 <- read.table(file="SRR306539_neocortex2.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
cortex2$sample <- "cortex2"
cortex3 <- read.table(file="SRR306540_neocortex3.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
cortex3$sample <- "cortex3"

cerebellum1 <- read.table(file="SRR306526_cerebellum1.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
cerebellum1$sample <- "cerebellum1"
cerebellum2 <- read.table(file="SRR306527_cerebellum2.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
cerebellum2$sample <- "cerebellum2"
cerebellum3 <- read.table(file="SRR306528_cerebellum3.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
cerebellum3$sample <- "cerebellum3"

purkinje1 <- read.table(file="SRR306529_purkinje1.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
purkinje1$sample <- "purkinje1"
purkinje2 <- read.table(file="SRR306530_purkinje2.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
purkinje2$sample <- "purkinje2"
purkinje3 <- read.table(file="SRR306531_purkinje3.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
purkinje3$sample <- "purkinje3"

camk2a1 <- read.table(file="SRR306532_camk2a1.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
camk2a1$sample <- "camk2a1"
camk2a2 <- read.table(file="SRR306533_camk2a2.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
camk2a2$sample <- "camk2a2"
camk2a3 <- read.table(file="SRR306534_camk2a1.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
camk2a3$sample <- "camk2a3"

gad21 <- read.table(file="SRR306535_gad21.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
gad21$sample <- "gad21"
gad22 <- read.table(file="SRR306536_gad22.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
gad22$sample <- "gad22"
gad23 <- read.table(file="SRR306537_gad23.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
gad23$sample <- "gad23"

PV1 <- read.table(file="SRR306541_PV1.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
PV1$sample <- "PV1"
PV2 <- read.table(file="SRR306542_PV1.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
PV2$sample <- "PV2"

SST1 <- read.table(file="SRR306543_SST1.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
SST1$sample <- "SST1"
SST2 <- read.table(file="SRR306544_SST2.fastq.collapsed", header = FALSE, fill = FALSE, sep =  "")
SST2$sample <- "SST2"

brain <- rbind(cortex1,cortex2,cortex3,cerebellum1,cerebellum2,cerebellum3)
rm(cortex1,cortex2,cortex3,cerebellum1,cerebellum2,cerebellum3)
brain <- rbind(brain,camk2a1,camk2a2,camk2a3, gad21, gad22, gad23, purkinje1,purkinje2,purkinje3)
rm(camk2a1,camk2a2,camk2a3, gad21, gad22, gad23, purkinje1,purkinje2,purkinje3)
brain <- rbind(brain,SST1,SST2,PV1,PV2)
rm(SST1,SST2,PV1,PV2)

motif_mir9_5p <- "TGGTTATCTAGCT"
brain$mir95p <- regexpr(motif_mir9_5p, brain$V1)

brain$counts <- brain$V1
brain$counts <- gsub(pattern = "A",replacement = "",x = brain$counts)
brain$counts <- gsub(pattern = "T",replacement = "",x = brain$counts)
brain$counts <- gsub(pattern = "G",replacement = "",x = brain$counts)
brain$counts <- gsub(pattern = "C",replacement = "",x = brain$counts)
brain$counts <- gsub(pattern = "N",replacement = "",x = brain$counts)
brain$counts <- gsub(pattern = " ",replacement = "",x = brain$counts)
brain$counts <- as.numeric(brain$counts)

mir95p <- brain[which(brain$mir95p!="-1"),]

mir95p_agg <- aggregate(mir95p$counts, by=list(SAMPLES=mir95p$sample, START=mir95p$mir95p), FUN=sum)
mir95p_total <- aggregate(mir95p$counts, by=list(SAMPLES=mir95p$sample), FUN=sum)
mir95p_agg <- merge(mir95p_agg, mir95p_total, by="SAMPLES")
mir95p_agg$percent <- mir95p_agg$x.x/mir95p_agg$x.y*100

reads_total <- aggregate(brain$counts, by=list(SAMPLES=brain$sample), FUN=sum)
reads_total <- merge(reads_total, mir95p_total, by="SAMPLES")
reads_total$CPM <- reads_total$x.y/reads_total$x.x*1000000


motif_mir9_3p <- "AGCTAGATAACCG"
brain$mir93p <- regexpr(motif_mir9_3p, brain$V2)

brain$counts <- brain$V1
brain$counts <- as.numeric(brain$counts)

mir93p <- brain[which(brain$mir93p!="-1"),]

mir93p_agg <- aggregate(mir93p$counts, by=list(SAMPLES=mir93p$sample, START=mir93p$mir93p), FUN=sum)
mir93p_total <- aggregate(mir93p$counts, by=list(SAMPLES=mir93p$sample), FUN=sum)
mir93p_agg <- merge(mir93p_agg, mir93p_total, by="SAMPLES")
mir93p_agg$percent <- mir93p_agg$x.x/mir93p_agg$x.y*100

reads_total <- aggregate(brain$counts, by=list(SAMPLES=brain$sample), FUN=sum)
reads_total <- merge(reads_total, mir93p_total, by="SAMPLES")
reads_total$CPM <- reads_total$x.y/reads_total$x.x*1000000





