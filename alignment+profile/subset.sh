
marine() {
	awk 'NR==FNR{_[$1]=1}NR>FNR&&_[$1]' - <(awk -F, '$3~/marine/{print$1}' ../datasets/habitats.txt) | xargs samtools faidx nbla.mafft-trimmed.fasta
}

samtools faidx nbla.mafft-trimmed.fasta
awk '/Viruses/{print$1}'               ../uniprot.tax.override   | marine > nbla.mafft-viruses.fasta
awk -F, '$3=="Cyanobacteria"{print$1}' ../datasets/taxonomy1.txt | marine > nbla.mafft-cyanobacteria.fasta
