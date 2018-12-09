library(ggplot2)

hmmsearch <- "NF90-pfam.out"
nbla  <- read.table(pipe(paste("awk", "'{for(i=24;i<=NF;i++)$23=$23\" \"$i;NF=23}1'", "OFS=\\\\t", hmmsearch)), comment.char = "#", sep = "\t", quote = "")
interpro <- read.table("NF90-interpro-nbla.tab", sep="\t", quote="", header = T)

cols <- c('target.name', 'accession', 'tlen', 'query.name', 'accession2', 'qlen', 'full.e.value', 'full.score', 'full.bias', 'domain.num', 'domain.of', 'domain.c.evalue', 'domain.i.value', 'domain.score', 'domain.bias', 'hmm.from', 'hmm.to', 'ali.from', 'ali.to', 'env.from', 'env.to', 'acc', 'description.of.target')

names(nbla)[1:length(cols)] <- cols

nbla <- subset(nbla, !grepl("UniRef90_UPI", target.name))
nbla <- nbla[order(-nbla$domain.score),]
nbla <- nbla[!duplicated(nbla$target.name),]

nbla$ali.len <- nbla$ali.to - nbla$ali.from + 1
nbla$hmm.len <- nbla$hmm.to - nbla$hmm.from + 1
nbla$env.len <- nbla$env.to - nbla$env.from + 1

nbla <- merge(nbla, interpro, by.x = "target.name", by.y = "Cluster.ID", all.x = T)

nbla$interpro <- !is.na(nbla$Cluster.name)

nbla$Inclusion <- with(nbla, ifelse(hmm.len > 44, "Profile matches", "Rejected sequences"))
nbla$Category <- with(nbla, ifelse(interpro, 1, ifelse(grepl("gene", target.name), 2, 3)))
nbla$Category <- as.factor(nbla$Category)
levels(nbla$Category) <- c("Annotated in InterPro", "Phage", "Other sequences")

ggsave(paste0(hmmsearch, "-hmmscan.svg"),
	ggplot(nbla, aes(x = hmm.len, y = log10(domain.i.value), color = Category)) + geom_point() +
	xlab("Alignment length") + ylab("log e-value") + theme_bw()
)

