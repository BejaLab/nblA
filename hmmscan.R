library(ggplot2)
library(ggrepel)

nbla  <- read.table("NF50-pfam.out",  skip = 3)
phage <- read.table("phage-pfam.out", skip = 3)
interpro <- read.table("NF50-interpro-nbla.tab", sep="\t", quote="", header=T)
final    <- read.table("NF50-pfam+interpro.tab", sep="\t", quote="", header=T)

cols <- c('target.name','accession','tlen','query.name','accession2','qlen','full.e.value','full.score','full.bias','domain.num','domain.of','domain.c.evalue','domain.i.value','domain.score','domain.bias','hmm.from','hmm.to','ali.from','ali.to','env.from','env.to','acc')

nbla$phage  <- F
phage$phage <- T
nbla <- rbind(nbla, phage)
colnames(nbla)[1:length(cols)] <- cols

nbla <- nbla[order(nbla$domain.i.value),]
nbla <- nbla[!duplicated(nbla$query.name),]
nbla$ali.len <- nbla$ali.to - nbla$ali.from + 1

nbla <- merge(nbla, interpro, by.x = "query.name", by.y = "Cluster.ID", all.x = T)
nbla <- merge(nbla, final,    by.x = "query.name", by.y = "Cluster.ID", all.x = T)

nbla$interpro <- !is.na(nbla$Cluster.name.x)
nbla$final    <- !is.na(nbla$Cluster.name.y)

nbla$Inclusion <- with(nbla, ifelse(domain.i.value < 1e-8 | domain.i.value < 1e-5 & ali.len > 42, "Profile matches", ifelse(final, "Additional blast matches", "Rejected sequences")))
nbla$Inclusion <- factor(nbla$Inclusion, levels = c("Rejected sequences", "Profile matches", "Additional blast matches"))
#nbla$label <- as.character(nbla$query.name)
nbla$Category <- with(nbla, ifelse(phage, "Phage MAGs", ifelse(interpro, "Annotated in InterPro", "Other sequences")))
nbla$Category <- factor(nbla$Category, levels = c("Annotated in InterPro", "Phage MAGs", "Other sequences"))

#nbla[!(nbla$query.name %in% c("gene1", "gene1b")),"label"] <- ""

ggsave("hmmscan-results.svg",
	ggplot(nbla, aes(x = ali.len, y = log10(domain.i.value), color = Category, shape = Inclusion)) +
	geom_point() +
	#geom_label_repel(aes(label = label)) +
	xlab("Alignment length") + ylab("log e-value") +
	theme_bw()
)
