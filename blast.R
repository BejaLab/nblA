library(ggplot2)
library(ggrepel)

pepe  <- read.table("pepe-uniref50.tab", sep="\t", quote="", header=T)
blast <- read.table("NF50-pfam+interpro.blast", col.names = c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'))
queries <- levels(blast$qseqid)
blast <- subset(blast, !(sseqid %in% queries))
blast <- blast[rev(order(blast$bitscore)),]
blast <- blast[!duplicated(blast$sseqid),]

target <- c('UniRef50_Q05RZ2','UniRef50_Q3AZR7','UniRef50_A0A1Z9W551','UniRef50_Q064K8','UniRef50_A0A1Z8PBX3','UniRef50_Q0I8F0','UniRef50_A0A2D6Y599','UniRef50_A0A0H4BCM9','UniRef50_A0A164D6M4','UniRef50_A0A076H7M8','UniRef50_A0A2D5RAY9','UniRef50_A0A1Z9RDX9','UniRef50_G4FIY3','UniRef50_A0A076H0N1','UniRef50_A0A163W0S2','UniRef50_A0A076HGV9','UniRef50_A0A1Z9KU07','UniRef50_A0A2E0KFG2','UniRef50_A0A076HHW8','UniRef50_A0A1Z8W9F8','UniRef50_A0A163Y795','UniRef50_W0GUF6','UniRef50_A0A164B425','UniRef50_A0A2D6FH08','UniRef50_A0A2D5RB34','UniRef50_A0A081GI73','UniRef50_A0A2T1ET24','UniRef50_K9P6Y6','UniRef50_A0A316JT77','UniRef50_A0A2E8TDC2','UniRef50_A0A1Z9JBD0','UniRef50_A0A0H5PUN8','UniRef50_A0A2R7TRJ3','UniRef50_A4CTT0','UniRef50_A0A1Z9WA43','UniRef50_A5GKC6','UniRef50_A0A081GLQ4','UniRef50_A0A326QGI7','UniRef50_B5IM81','UniRef50_A5GQX1','UniRef50_A0A2D6Y678','UniRef50_A3YW23','UniRef50_A0A316JM73','UniRef50_A5GK08','UniRef50_Q05R15','UniRef50_A0A2G4HMK1','UniRef50_A0A182AP08','UniRef50_A0A2E0KI01','UniRef50_A0A2E0AGH7','UniRef50_A0A1Z8ND19','UniRef50_A4CXU5','UniRef50_A0A251Z8Z4','UniRef50_A0A076GWZ4','UniRef50_A0A182ATG2','UniRef50_B5IL66','UniRef50_A0A2P7EHR8','UniRef50_A0A1J0PC72','UniRef50_A0A2P7EBP5','UniRef50_A0A251Z8A2','UniRef50_A0A326QD80','UniRef50_A0A2D8TVK8','UniRef50_A0A2D8U0Y7')

blast$qlen <- blast$qend - blast$qstart + 1
blast$slen <- blast$send - blast$sstart + 1

blast <- merge(blast, pepe, by.x = "sseqid", by.y = "Cluster.ID", all.x = T)
blast$group <- ifelse(blast$sseqid %in% target, "Target", ifelse(is.na(blast$Cluster.name), "New", "Pepe"))

ggplot(blast, aes(x = qlen, y = log10(evalue), color = group)) + geom_point()
