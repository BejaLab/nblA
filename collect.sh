
uniprotList() {
	local query=$(sort | uniq | xargs)
	local data="from=$1&to=$2&format=$3&query=$query"
	curl -L -d "$data" -X POST https://www.uniprot.org/uploadlists/
}

uniprotTax() {
	local query=$(sort | uniq | xargs | sed -e 's/ /+OR+/g')
	local data="query=$query&columns=id,lineage(ALL)&format=tab"
	curl -L -d "$data" https://www.uniprot.org/uniprot/
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
[ -s uniprot.tab ]   || cut -f6 NF50-pfam+interpro.tab | tail -n+2 | tr \; \\n | tr -d ' ' | grep -v ^UPI | uniprotList ACC ACC tab > uniprot.tab
[ -s uniprot.fasta ] || cut -f2 uniprot.tab | tail -n+2 | uniprotList ACC ACC fasta > uniprot.fasta
[ -s uniprot.tax ]   || cut -f2 uniprot.tab | tail -n+2 | uniprotTax > uniprot.tax

cut -f2 uniprot.tax | tail -n+2 | awk -F, '{t=$3}$3~/group/{t=$4}$4~/group/{t=$5}$1~/Virus/{t=$1} {print t}'
cut -f2 uniprot.tax | tail -n+2 | sed -e 's/(.*)//' -e 's/, /,/g' | awk -F, '{t1=$3;t2=$4;t3=""}$3~/group/{t1=$4;t2=$5}$4~/group/{t1=$5;t2=$6}$1~/Virus/{t1=$1;t2=$4}$5~/viridae/{t2=$4"-"$5} !t2||tolower(t1 t2)~/unclassified|uncultured|candidat/{t2=""} {print t1,t2,t3}'

cat <<- EOF > taxonomy1.dataset
	DATASET_TEXT
	SEPARATOR COMMA
	DATASET_LABEL,Taxonomy - level 2
	DATA
	$(tail -n+2 uniprot.tax| tr \\t , | sed -e 's/(.*)//' -e 's/, /,/g' | awk -F, '{t=$4}$4~/group/{t=$5}$5~/group/{t=$6}$2~/Virus/{t=$2} {print $1,t}' OFS=,)
EOF

cat <<- EOF > taxonomy2.dataset
	DATASET_TEXT
	SEPARATOR COMMA
	DATASET_LABEL,Taxonomy - level 2
	DATA
	$(tail -n+2 uniprot.tax| tr \\t , | sed -e 's/(.*)//' -e 's/, /,/g' | awk -F, '{t=$5}$4~/group/{t=$6}$5~/group/{t=$7}$1~/Virus/{t=$5}$6~/viridae/{t=$5"-"$6} !t||tolower(t)~/unclassified|uncultured|candidat/{t=""} {print $1,t}' OFS=,)
EOF
