
uniprotList() {
	local query=$(sort | uniq | xargs)
	local data="from=$1&to=$2&format=$3&query=$query"
	[ -z "$4" ] || data="$data&columns=$4"
	curl -L -d "$data" -X POST https://www.uniprot.org/uploadlists/
}

emblENA() {
	curl -L "https://www.ebi.ac.uk/ena/data/view/$1%26display%3Dxml"
}

biosamples() {
	local db=biosample
	xargs -n20 | sed -e 's/ /[accession]+OR+/g' -e s/+OR+$// | while read query; do
		query=$(curl -d "db=$db&term=$query" https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi | xmlstarlet fo -D | xmlstarlet sel -T -t -m /eSearchResult/IdList/Id -v . -n | awk 1 ORS=, | sed s/,$//)
		curl -d "db=$db&id=$query" https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi | xmlstarlet fo -D | xmlstarlet sel -T -t -m /BioSampleSet/BioSample -v @accession -o $'\t' \
		-v "Attributes/Attribute[@attribute_name='env_biome']"        -o / \
		-v "Attributes/Attribute[@attribute_name='env_feature']"      -o / \
		-v "Attributes/Attribute[@attribute_name='env_material']"     -o / \
		-v "Attributes/Attribute[@attribute_name='geo_loc_name']"     -o / \
		-v "Attributes/Attribute[@attribute_name='isolation_source']" -o / \
		-v "Attributes/Attribute[@attribute_name='lat_lon']" -n
	done
}

parseENA() {
	sed -E 's/^[acgtn ]+$//' | \
		xmlstarlet sel -T -t -m /ROOT/entry --var lb -n --break -v @accession -o $'\t' -v 'xref[@db="BioSample"]/@id' -o $'\t' -v 'translate(.//qualifier[@name="isolation_source"]/value,$lb,"")' -n
}

[ -s PF04485.hmm ] || curl http://pfam.xfam.org/family/PF04485/hmm > PF04485.hmm && hmmpress PF04485.hmm
[ -s NF90-pfam.out ]     || cat phage.fasta /opt/uniprot/uniref90.fasta | hmmsearch       --cpu 4 --domtblout NF90-pfam.out     PF04485.hmm - &> NF90-pfam.log
[ -s NF90-pfam-T14.out ] || cat phage.fasta /opt/uniprot/uniref90.fasta | hmmsearch -T 14 --cpu 4 --domtblout NF90-pfam-T14.out PF04485.hmm - &> NF90-pfam-T14.log

[ -s interpro-nbla.fasta ] || curl http://www.ebi.ac.uk/interpro/entry/IPR036904/proteins-matched?export=fasta > interpro-nbla.fasta
[ -s NF90-interpro-nbla.tab ] || grep '>' interpro-nbla.fasta | tr -d '>' | uniprotList ACC NF90 tab > NF90-interpro-nbla.tab

[ -s NF90-pfam-filtered.tab ] || awk '$17-$16>43{print$1}' NF90-pfam.out | grep UniRef90 | sort | uniq | uniprotList NF90 NF90 tab > NF90-pfam-filtered.tab
[ -s NF90-pfam-NF50.tab ] || cut -f2 NF90-pfam-filtered.tab NF90-interpro-nbla.tab | cut -d_ -f2 | grep -v Cluster | uniprotList ACC NF50 tab id,members > NF90-pfam-NF50.tab
[ -s uniprot-pfam.fasta ] || csvtool -t TAB namedcol 'Cluster members' NF90-pfam-NF50.tab | tr \; \\n | tr -d ' ' | grep -v Cluster | uniprotList ACC ACC fasta > uniprot-pfam.fasta

#MAFFT_BINARIES= mafft <(cat uniprot-pfam.fasta phage.fasta) > nbla1.mafft
#MAFFT_BINARIES= mafft uniprot-pfam.fasta > nbla2.mafft
#[ -s nbla1.hmm ] || hmmbuild nbla1.hmm nbla1.mafft && hmmpress nbla1.hmm
#[ -s nbla2.hmm ] || hmmbuild nbla2.hmm nbla2.mafft && hmmpress nbla2.hmm
#[ -s NF90-nbla1.out ] || cat phage.fasta /opt/uniprot/uniref90.fasta | hmmsearch --cpu 4 --domtblout NF90-nbla1.out nbla1.hmm - &> nbla1.log
#[ -s NF90-nbla2.out ] || cat phage.fasta /opt/uniprot/uniref90.fasta | hmmsearch --cpu 4 --domtblout NF90-nbla2.out nbla2.hmm - &> nbla2.log

samtools faidx uniprot-pfam.fasta
[ -s uniprot-pfam.blast ] || \
	parallel --colsep \\t -j 12 "samtools faidx uniprot-pfam.fasta {1} | blastp -query /dev/stdin -outfmt 6 -db /opt/uniprot/uniref90/uniref90" < uniprot-pfam.fasta.fai > uniprot-pfam.blast

samtools faidx phage.fasta
[ -s phage.blast ] || parallel --colsep \\t -j 12 "samtools faidx phage.fasta {1} | blastp -query /dev/stdin -outfmt 6 -db /opt/uniprot/uniref90/uniref90" < phage.fasta.fai > phage.blast

[ -s NF90-pfam.tab ] || awk '$11<0.03{print$2}' uniprot-pfam.blast phage.blast | sort | uniq | uniprotList NF90 NF90 tab > NF90-pfam.tab

[ -s NF90-all-NF50.tab ] || cut -f2 NF90-pfam.tab | cut -d_ -f2 | grep -v Cluster | uniprotList ACC NF50 tab id,members > NF90-all-NF50.tab

[ -s uniprot.tab ]   || csvtool -t TAB namedcol 'Cluster members' NF90-all-NF50.tab | tail -n+2 | tr \; \\n | tr -d ' ' | grep -v ^UPI | \
	uniprotList ACC ACC tab "id,organism-id,database(EMBL),proteome,lineage(ALL)" > uniprot.tab
[ -s uniprot.fasta ] || csvtool -t TAB namedcol Entry uniprot.tab | tail -n+2 | uniprotList ACC ACC fasta > uniprot.fasta

#cut -f1 uniprot.tab | tail -n+2 | uniprotList ACC EMBL_ID tab > uniprot-embl.tab

#if [ ! -s uniprot-embl.source ]; then
#	printf "%s\t%s\t%s\n" accession biosample isolation_source > uniprot-embl.source
#	cut -f2 uniprot-embl.tab | tail -n+2 | xargs -n 10 | tr ' ' , | while read query; do emblENA "$query" | parseENA; done >> uniprot-embl.source
#fi

#comm -23 <(cut -f2 uniprot-embl.tab | tail -n+2 | sort) <(cut -f1 uniprot-embl.source | tail -n+2 | sort) | while read query; do emblENA "$query" | parseENA; done >> uniprot-embl.source

#cut -f2 uniprot-embl.source | tail -n+2 | biosamples > uniprot-embl.biosample

# awk -F\\t 'FNR>1&&FILENAME=="uniprot.tax"{e=em[$1];print $1,$2,sc[e],sam[e],env[sam[e]]}FILENAME=="uniprot-embl.tab"{em[$1]=$2}FILENAME=="uniprot-embl.source"{sam[$1]=$2;sc[$1]=$3}FILENAME=="uniprot-embl.biosample"{env[$1]=$2}' OFS=\\t uniprot-embl.tab uniprot-embl.source uniprot-embl.biosample uniprot.tax > habitats.tsv
