library(ggplot2)
library(ggrepel)

nbla  <- read.table("NF50-pfam.out",  skip = 3)
phage <- read.table("phage-pfam.out", skip = 3)
interpro <- read.table("interpro-nbla-NF50.tab", sep="\t", quote="", header=T)

cols <- c('target.name','accession','tlen','query.name','accession2','qlen','full.e.value','full.score','full.bias','domain.num','domain.of','domain.c.evalue','domain.i.value','domain.score','domain.bias','hmm.from','hmm.to','ali.from','ali.to','env.from','env.to','acc')

nbla$phage  <- F
phage$phage <- T
nbla <- rbind(nbla, phage)
colnames(nbla)[1:length(cols)] <- cols

nbla$ali.len <- nbla$ali.to - nbla$ali.from + 1
nbla <- merge(nbla, interpro, by.x = "query.name", by.y = "Cluster.ID", all.x = T)

nbla$label <- as.character(nbla$query.name)
nbla$group <- ifelse(nbla$phage, "Phage", ifelse(is.na(nbla$Cluster.name), "New", "Interpro"))
nbla[!nbla$phage & !(nbla$query.name %in% c('UniRef50_A0A316JT77','UniRef50_A0A164B425')),"label"] <- ""

ggplot(nbla, aes(x = ali.len, y = log10(domain.i.value), color = group)) + geom_point() + geom_label_repel(aes(label = label))
