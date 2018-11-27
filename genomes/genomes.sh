
# ncbi query: ("complete genome"[Assembly Level]) AND cyanobacteria[organism]

cat fasta/*_genomic.fna | getorf -find 0 -circular -minsize 100 -table 11 -filter | cat - fasta/*_protein.faa > target.fa
find fasta/*.f*a -exec samtools faidx {} \;
awk '{print $1,FILENAME}' fasta/*.fai > index.txt

makeblastdb -dbtype prot -in target.fa -out target

hmmsearch --cpu 4 --domtblout target.out ../PF04485.hmm target.fa &> target.log

parallel --colsep \\t -j 12 "samtools faidx ../uniprot-pfam.fasta {1} | blastp -query /dev/stdin -outfmt 6 -db target" < ../uniprot-pfam.fasta.fai > target.blast
