
source params.cfg

mkdir hmm uniprot

uniprotList() {
	local query=$(sort | uniq | xargs)
	local data="from=$1&to=$2&format=$3&query=$query"
	[ -z "$4" ] || data="$data&columns=$4"
	curl -L -d "$data" -X POST https://www.uniprot.org/uploadlists/
}

# Get the nblA hmm profile
curl http://pfam.xfam.org/family/PF04485/hmm > hmm/PF04485.hmm && hmmpress hmm/PF04485.hmm

# Search uniref90 for PF04485 matches
cat phage.fasta "$UNIREF90/uniref90.fasta" | hmmsearch --cpu 4 --domtblout uniprot/01_NF90-pfam.out PF04485.hmm - &> uniprot/01_NF90-pfam.log
awk '$17-$16>43{print$1}' uniprot/01_NF90-pfam.out | grep UniRef90 | sort | uniq | uniprotList NF90 NF90 tab > uniprot/01_NF90-pfam-filtered.tab

# Download proteins matched to IPR036904 and their Uniref90 cluster proteins
curl http://www.ebi.ac.uk/interpro/entry/IPR036904/proteins-matched?export=fasta > uniprot/02_interpro-nbla.fasta
seqkit seq -n uniprot/02_interpro-nbla.fasta | uniprotList ACC NF90 tab > uniprot/02_NF90-interpro-nbla.tab

# Expand the space by matching to Uniref50 
cut -f2 uniprot/01_NF90-pfam-filtered.tab uniprot/02_NF90-interpro-nbla.tab | cut -d_ -f2 | grep -v Cluster | uniprotList ACC NF50 tab id,members > uniprot/03_NF90-pfam-NF50.tab
csvtool -t TAB namedcol 'Cluster members' uniprot/03_NF90-pfam-NF50.tab | tr \; \\n | tr -d ' ' | grep -v Cluster | uniprotList ACC ACC fasta > uniprot/03_uniprot-pfam.fasta

# Blast the proteins agains Uniref90
blastp -num_threads "$CPUs" -query uniprot/03_uniprot-pfam.fasta -outfmt 6 -db "$UNIREF90" -out uniprot/04_uniprot-pfam.blast
blastp -num_threads "$CPUs" -query phage.fasta -outfmt 6 -db "$UNIREF90/uniref90" -out uniprot/04_phage.blast

# Filter the blast results by E-value threshold of 0.03
awk '$11<0.03{print$2}' uniprot/04_{uniprot-pfam,phage}.blast | sort | uniq | uniprotList NF90 NF90 tab > uniprot/05_NF90-pfam.tab

# Expand the matches again via Uniref50
# The final fasta file is uniprot.fasta
cut -f2 uniprot/05_NF90-pfam.tab | cut -d_ -f2 | grep -v Cluster | uniprotList ACC NF50 tab id,members > uniprot/06_NF90-all-NF50.tab
csvtool -t TAB namedcol 'Cluster members' uniprot/06_NF90-all-NF50.tab | tail -n+2 | tr \; \\n | tr -d ' ' | grep -v ^UPI | \
	uniprotList ACC ACC tab "id,organism-id,database(EMBL),proteome,lineage(ALL)" > uniprot/07_uniprot.tab
csvtool -t TAB namedcol Entry uniprot/07_uniprot.tab | tail -n+2 | uniprotList ACC ACC fasta > uniprot.fasta
