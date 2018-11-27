library(ggplot2)
library(stringr)

cmd  <- "sed -E 's/: +/\\t/' fasta/*_assembly_report.txt | tr -d \\\\r | awk -F\\\\t '/^#/{_[$1]=$2;next}{print$0,_[\"# Organism name\"],_[\"# Taxid\"],_[\"# GenBank assembly accession\"]}' OFS=\\\\t"
cols <- c("Sequence.Name", "Sequence.Role", "Assigned.Molecule", "Assigned.Molecule.Location.Type", "GenBank.Accn", "Relationship", "RefSeq.Accn", "Assembly.Unit", "Sequence.Length", "UCSC.style.name", "Organism.Name", "Taxid", "assembly")
genomic <- read.table(pipe(cmd), sep = "\t", quote = "", col.names = cols, stringsAsFactors = F)
assemblies <- subset(genomic, !duplicated(assembly), c("Organism.Name", "Taxid", "assembly"))

cmd  <- "cat fasta/*_feature_table.txt"
cols <- c("feature", "class", "assembly", "assembly_unit", "seq_type", "chromosome", "genomic_accession", "start", "end", "strand", "product_accession", "non-redundant_refseq", "related_accession", "name", "symbol", "GeneID", "locus_tag", "feature_interval_length", "product_length", "attributes")
features <- read.table(pipe(cmd), sep = "\t", quote = "", col.names = cols, stringsAsFactors = F)

cmd  <- paste("awk", "'{for(i=24;i<=NF;i++)$23=$23\" \"$i;NF=23}1'", "OFS=\\\\t", "target.out")
cols <- c('target.name', 'accession', 'tlen', 'query.name', 'accession2', 'qlen', 'full.e.value', 'full.score', 'full.bias', 'domain.num', 'domain.of', 'domain.c.evalue', 'domain.i.value', 'domain.score', 'domain.bias', 'hmm.from', 'hmm.to', 'ali.from', 'ali.to', 'env.from', 'env.to', 'acc', 'description.of.target')
out <- read.table(pipe(cmd), comment.char = "#", sep = "\t", quote = "", col.names = cols, stringsAsFactors = F)

out <- out[order(-out$domain.score),]
out <- out[!duplicated(out$target.name),]

out$hmm.len <- with(out, hmm.to - hmm.from + 1)
out <- subset(out, hmm.len > 43)

cols <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')
blast <- read.table("target.blast", col.names = cols, stringsAsFactors = F)
blast <- blast[order(-blast$bitscore),]
blast <- blast[!duplicated(blast$sseqid),]
blast <- subset(blast, evalue < 0.03)

matches <- c(out$target.name, blast$sseqid)
matches <- matches[!duplicated(matches)]
matches <- str_remove(matches, "_\\d+$")
matches <- data.frame(table(matches))

proteins  <- merge(matches, features, all.x = T, by.x = "matches", by.y = "product_accession")
proteins <- aggregate(Freq ~ assembly, proteins, sum)
proteins <- merge(proteins, assemblies, all = T)
proteins$type <- 2

orfs      <- merge(matches, genomic, by.x = "matches", by.y = "GenBank.Accn")
orfs      <- aggregate(Freq ~ assembly, orfs, sum)
orfs      <- merge(orfs, assemblies, all = T)
orfs$type <- 3

nblas <- subset(features, class == "with_protein" & (grepl("phycobilisome degrdation protein|nbla", tolower(name)) | grepl("nbla", tolower(symbol))))
nblas <- aggregate(product_accession ~ assembly, nblas, length)
nblas$Freq <- nblas$product_accession
nblas <- merge(nblas, assemblies, all = T)
nblas$type <- 1

cols <- c("assembly","Freq","type")
freqs <- rbind(proteins[,cols], orfs[,cols], nblas[,cols])
genomes <- merge(assemblies, freqs, all = T)
genomes[is.na(genomes)] <- 0
genomes$type <- as.factor(genomes$type)
levels(genomes$type) <- c("Annotated as nblA", "Predicted proteins", "All ORFs")

prochlo <- subset(genomes,  grepl("Prochlorococcus", Organism.Name))
genomes <- subset(genomes, !grepl("Prochlorococcus", Organism.Name))

ggsave("genomes.svg",
	ggplot(genomes, aes(x = Freq, fill = type)) +
		geom_bar(position = position_dodge2(preserve = "single")) +
		xlab("Number of nblA genes") + ylab("Number of genomes") +
		scale_x_continuous(breaks = 0:7) +
		scale_fill_manual(values = c("#cccccc","#969696","#525252")) +
		theme_bw() + theme(legend.position = c(0.8,0.8), legend.title = element_blank()),
	width = 7, height = 5
)
