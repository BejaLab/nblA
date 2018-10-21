
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

[ -s phage-pfam.out ] || hmmscan --cpu 2 --domtblout phage-pfam.out NblA.hmm phage.fasta &> phage-pfam.log
[ -s NF50-pfam.out  ] || hmmscan --cpu 2 --domtblout NF50-pfam.out NblA.hmm /opt/uniprot/uniref50.fasta &> NF50-pfam.log

[ -s interpro-nbla.fasta ] || curl http://www.ebi.ac.uk/interpro/entry/IPR036904/proteins-matched?export=fasta > interpro-nbla.fasta
[ -s NF50-interpro-nbla.tab ] || grep '>' interpro-nbla.fasta | tr -d '>' | uniprotList ACC NF50 tab > NF50-interpro-nbla.fasta

[ -s NF50-pfam-filtered.tab ] || awk '$13<1e-8||$13<1e-5&&$19-$18>42{print$4}' NF50-pfam.out | grep UniRef50 | uniprotList NF50 NF50 tab > NF50-pfam-filtered.tab
[ -s NF50-pfam+interpro.fasta ] || cut -f2 NF50-pfam-filtered.tab NF50-interpro-nbla.fasta | grep -v ^Cluster | sort | uniq | uniprotList NF50 NF50 fasta > NF50-pfam+interpro.fasta

samtools faidx NF50-pfam+interpro.fasta
[ -s NF50-pfam+interpro.blast ] || parallel --colsep \\t -j 12 "samtools faidx NF50-pfam+interpro.fasta {1} | blastp -query /dev/stdin -outfmt 6 -db /opt/uniprot/uniref50/uniref50" < NF50-pfam+interpro.fasta.fai > NF50-pfam+interpro.blast

[ -s NF50-pfam+interpro.tab ] || awk '$11<0.03{print$1;print$2}' NF50-pfam+interpro.blast | sort | uniq | uniprotList NF50 NF50 tab > NF50-pfam+interpro.tab
[ -s uniprot.tab ]   || cut -f6 NF50-pfam+interpro.tab | tail -n+2 | tr \; \\n | tr -d ' ' | grep -v ^UPI | uniprotList ACC ACC tab "id,lineage(ALL)" > uniprot.tab
[ -s uniprot.fasta ] || cut -f1 uniprot.tab | tail -n+2 | uniprotList ACC ACC fasta > uniprot.fasta

cut -f1 uniprot.tab | tail -n+2 | uniprotList ACC EMBL_ID tab > uniprot-embl.tab

if [ ! -s uniprot-embl.source ]; then
	printf "%s\t%s\t%s\n" accession biosample isolation_source > uniprot-embl.source
	cut -f2 uniprot-embl.tab | tail -n+2 | xargs -n 10 | tr ' ' , | while read query; do emblENA "$query" | parseENA; done >> uniprot-embl.source
fi

comm -23 <(cut -f2 uniprot-embl.tab | tail -n+2 | sort) <(cut -f1 uniprot-embl.source | tail -n+2 | sort) | while read query; do emblENA "$query" | parseENA; done >> uniprot-embl.source

cut -f2 uniprot-embl.source | tail -n+2 | biosamples > uniprot-embl.biosample

awk -F\\t 'FNR>1&&FILENAME=="uniprot.tax"{e=em[$1];print $1,$2,sc[e],sam[e],env[sam[e]]}FILENAME=="uniprot-embl.tab"{em[$1]=$2}FILENAME=="uniprot-embl.source"{sam[$1]=$2;sc[$1]=$3}FILENAME=="uniprot-embl.biosample"{env[$1]=$2}' OFS=\\t uniprot-embl.tab uniprot-embl.source uniprot-embl.biosample uniprot.tax > habitats.tsv
