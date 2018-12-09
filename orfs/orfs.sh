
proteomes() {
	local query=$1
	curl -s "https://www.uniprot.org/proteomes/?query=($query)+redundant:no&format=tab&force=true&columns=id,name,organism-id,lineage,proteincount,proteome-components,assembly"
}

ena() {
	local id=$1
	local format=$2
	curl -sL "https://www.ebi.ac.uk/ena/data/view/$id&display=$format"
}

ena_sequence() {
	local id=$1
	local format=$2
	curl -sL "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ena_sequence&id=$id&format=$format&style=raw"
}

url_link() {
	local xml=$1
	local label=$2
	xmlstarlet sel -t -m "//URL_LINK[LABEL='$label']" -v URL "$xml"
}

tsvcols() {
	local file=$1
	local cols=$2
	csvtool -t TAB -u TAB namedcol "$cols" "$file" | tail -n+2
}

nbla_list() {
	cat <(awk '$17-$16>43{print$1}' target.out) <(awk '$11<0.0001{print$2}' target.blast) | sort | uniq
}

if [ ! -s proteomes.tab ]; then
	#proteomes taxonomy:1117           > proteomes.tab # cyanobacteria
	proteomes '((nostoc+OR+synechococcus+OR+plaktothrix+OR+microcystis+OR+phormidium+OR+anacystis+OR+lyngbya+OR+plectonema+OR+synechocystis)+(phage+OR+virus))+OR+cyanophage+OR+vB_NpeS-2AV2+OR+MedDCM-OCT-S09-C28' > proteomes.tab # cyanophages
fi

mkdir -p xml dat fa
tsvcols proteomes.tab 'Proteome ID,Genome assembly ID' | while read proteome assembly; do
	echo $proteome
	xml=xml/$proteome.xml
	dat=dat/$proteome.dat
	[ -s "$xml" ] || ena "$assembly" xml > "$xml"
	if [ ! -s "$dat" ]; then
		report=$(url_link "$xml" 'Sequence Report')
		if [ ! -z "$report" ]; then
			curl -sL "$report" | tsvcols - accession | while read accession; do
				ena_sequence "$accession" embl
			done > $dat
		else
			flat=$(url_link "$xml" 'WGS_SET_FLATFILE')
			if [ ! -z "$flat" ]; then
				curl -sL "$flat" | gzip -cd > "$dat"
			else
				xmlstarlet sel -t -m //CHROMOSOMES/CHROMOSOME -v @accession -n "$xml" | while read chr; do
					contigs=$(ena "$chr" xml | xmlstarlet sel -t -m //contig/range -v @accession -n)
					[ -z "$contigs" ] && ena "$chr" text
					for contig in $contigs; do
						ena_sequence "$contig" embl
					done
				done > $dat
			fi
		fi
		sleep 1
	fi
	if [ ! -s "$dat.fa" ] && grep ^CO $dat > /dev/null; then
		awk -F'[ ;]+' '$1=="AC"{print$2}' "$dat" | while read id; do
			ena_sequence "$id" fasta
		done > "$dat.fa"
	fi
	emboss_dat=$dat
	[ -s "$dat.fa" ] && emboss_dat=$dat.fa

	[ -s "fa/$proteome.fa" ] || perl merge.pl <(perl cds.pl "$dat") <(getorf -find 1 -table 11 -filter "$emboss_dat") 50 > "fa/$proteome.fa"
done

cat ../phage.fasta fa/*.fa > target.fa
[ -s target.out ] || hmmsearch --cpu 4 --domtblout target.out ../PF04485.hmm target.fa &> orfs-hmmsearch.log

[ -s target.pin ]   || makeblastdb -dbtype prot -in target.fa -out target
[ -s target.blast ] || parallel --colsep \\t -j 12 "samtools faidx ../uniprot-pfam.fasta {1} | blastp -query /dev/stdin -outfmt 6 -db target" < ../uniprot-pfam.fasta.fai > target.blast

samtools faidx target.fa
nbla_list | awk 'NR==FNR&&/>/{_[$1]=$0}NR>FNR&&_[">"$1]{print _[">"$1]}' target.fa -
nbla_list | awk 'NR==FNR{_[$1]=1}NR>FNR&&!_[$1]&&!/^gene/' ignored-orfs.txt - | awk -F\| 'NR==FNR&&/>/{_[$2]=1}NR>FNR&&!_[$0]' ../uniprot.fasta - | xargs samtools faidx target.fa > orf-nbla.fa
