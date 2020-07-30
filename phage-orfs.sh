
source params.cfg

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

mkdir -p orfs

# Get list of uniprot proteomes for cyanophages
[ -s orfs/proteomes.tab ] || proteomes '((nostoc+OR+synechococcus+OR+planktothrix+OR+microcystis+OR+phormidium+OR+anacystis+OR+lyngbya+OR+plectonema+OR+synechocystis+OR+prochlorococcus)+(phage+OR+virus))+OR+cyanophage+OR+vB_NpeS-2AV2' | awk -F\\t '!_[$3];{_[$3]=1}' > orfs/proteomes.tab

# Download nucleotide sequences
# and generate collections of ORFs - annotated and unannotated
mkdir -p orfs/xml orfs/dat orfs/fa
tsvcols orfs/proteomes.tab 'Proteome ID,Genome assembly ID' | while read proteome assembly; do
	echo $proteome
	xml=orfs/xml/$proteome.xml
	dat=orfs/dat/$proteome.dat
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

	[ -s "fa/$proteome.fa" ] || perl merge.pl <(perl cds.pl "$dat") <(getorf -find 1 -table 11 -filter "$emboss_dat") "$proteome" 50 > "fa/$proteome.fa"
done

cat phage.fasta orfs/fa/*.fa > orfs/target.fa
hmmsearch --cpu 4 --domtblout orfs/target.out hmm/PF04485.hmm orfs/target.fa &> orfs/orfs-hmmsearch.log

makeblastdb -dbtype prot -in orfs/target.fa
blastp -num_threads "$CPUs" -query uniprot/03_uniprot-pfam.fasta -outfmt 6 -db orfs/target.fa -out orfs/target.blast

samtools faidx orfs/target.fa
# nbla_list | awk 'NR==FNR&&/>/{_[$1]=$0}NR>FNR&&_[">"$1]{print _[">"$1]}' orfs/target.fa -
nbla_list | grep -v ^gene | awk -F\| 'NR==FNR&&/>/{_[$2]=1}NR>FNR&&!_[$0]' uniprot.fasta - | xargs samtools faidx orfs/target.fa > orfs/orf-nbla.fa
