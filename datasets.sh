
mkdir -p datasets

taxonomy() {
	csvtool -t TAB -u TAB namedcol "Entry,Taxonomic lineage (ALL)" uniprot/uniprot.tab | tail -n+2 | cat uniprot/uniprot.tax.override - | tr \\t , | sed -e 's/, /,/g' | grep -vf outliers.txt
}

cat <<- EOF > datasets/taxonomy1.txt
	DATASET_COLORSTRIP
	SEPARATOR COMMA
	DATASET_LABEL,Taxonomy
	COLOR_BRANCHES,0
	LEGEND_TITLE,Taxonomic groups
	LEGEND_SHAPES,$(awk -F, '{print 1}' ORS=, datasets/taxonomy1.list | sed s/,$//)
	LEGEND_LABELS,$(awk -F, '{print$2}' ORS=, datasets/taxonomy1.list | sed s/,$//)
	LEGEND_COLORS,$(awk -F, '{print$3}' ORS=, datasets/taxonomy1.list | sed s/,$//)
	DATA
	$(taxonomy | sed -e 's/(.*)//' | awk -F, 'NR==FNR{n[$1]=$2;c[$1]=$3;next} _[$1]{next}{_[$1]=1} { for (i=2;i<NF;i++) if (c[$i]) { print$1,c[$i],n[$i]; next } }' OFS=, datasets/taxonomy1.list -)
EOF

cat <<- EOF > datasets/taxonomy2.txt
	DATASET_COLORSTRIP
	SEPARATOR COMMA
	DATASET_LABEL,Taxonomy
	COLOR_BRANCHES,0
	LEGEND_TITLE,Taxonomic groups
	LEGEND_SHAPES,$(awk -F, '{print 1}' ORS=, datasets/taxonomy2.list | sed s/,$//)
	LEGEND_LABELS,$(awk -F, '{print$2}' ORS=, datasets/taxonomy2.list | sed s/,$//)
	LEGEND_COLORS,$(awk -F, '{print$3}' ORS=, datasets/taxonomy2.list | sed s/,$//)
	DATA
	$(taxonomy | sed -e 's/(.*)//' | \
		awk -F, '_[$1]{next}{_[$1]=1} {t=$4}$4~/group/{t=$5}$5~/group/{t=$6}$2~/Virus|Metagenomic/{t=$2} {print$1,t}' OFS=, | \
		awk -F, 'NR==FNR{c[$1]=$3;next} c[$2]{print $1,c[$2],$2}' OFS=, datasets/taxonomy2.list -)
EOF

ssconvert habitats.xls habitats.csv
cat <<- EOF > datasets/habitats.txt
	DATASET_COLORSTRIP
	SEPARATOR COMMA
	DATASET_LABEL,Habitats
	COLOR_BRANCHES,0
	LEGEND_TITLE,Habitats
	LEGEND_SHAPES,$(awk -F, '{print 1}' ORS=, datasets/habitats.list | sed s/,$//)
	LEGEND_LABELS,$(awk -F, '{print$2}' ORS=, datasets/habitats.list | sed s/,$//)
	LEGEND_COLORS,$(awk -F, '{print$3}' ORS=, datasets/habitats.list | sed s/,$//)
	DATA
	$(awk -F, 'NR==FNR{c[$1]=$3;next} c[$NF]{print $1,c[$NF],$NF}' OFS=, datasets/habitats.list habitats.csv | grep -vf outliers)
EOF

cat <<- EOF > datasets/organism.txt
	DATASET_TEXT
	SEPARATOR COMMA
	DATASET_LABEL,Organism
	DATA
	$(taxonomy | awk -F, '{print$1,$NF}' OFS=, | sed -e 's/(strain \(.*\))/\1/' -e 's/) (.*//g' -e 's/(.*//g' | awk -F, '!_[$1]{_[$1]=1; print$0,-1,"black","normal",1,0}' OFS=,)
EOF

cat <<- EOF > datasets/connections.txt
	DATASET_CONNECTION
	SEPARATOR COMMA
	DATASET_LABEL,Multiple NblA's
	DATA
	$(taxonomy | grep -v Rhodophyta | \
		awk -F, 'met[$1]{next}{met[$1]=1}!/sp.$/&&!/environmental/&&!/re-assigned/{_[$NF][$1]=1}END{for(sp in _) for (i1 in _[sp]) for (i2 in _[sp]) if (i1<i2) print i1,i2,1,"black",sp}' OFS=,)
EOF

cat <<- EOF > datasets/highlighted.txt
	DATASET_BINARY
	SEPARATOR COMMA
	DATASET_LABEL,Highlighted proteins
	FIELD_SHAPES,3
	FIELD_LABELS,star
	DATA
	$(awk '{print $1,1}' OFS=, datasets/highlighted.list)
EOF
