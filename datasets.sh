
mkdir -p datasets

cat <<- EOF > datasets/taxonomy.txt
	DATASET_COLORSTRIP
	SEPARATOR COMMA
	DATASET_LABEL,Taxonomy
	COLOR_BRANCHES,0
	LEGEND_TITLE,Taxonomic groups
	LEGEND_SHAPES,$(awk -F, '{print 1}' ORS=, datasets/taxonomy.list | sed s/,$//)
	LEGEND_LABELS,$(awk -F, '{print$2}' ORS=, datasets/taxonomy.list | sed s/,$//)
	LEGEND_COLORS,$(awk -F, '{print$3}' ORS=, datasets/taxonomy.list | sed s/,$//)
	DATA
	$(cut -f1,2 uniprot.tax | tail -n+2 | cat uniprot.tax.override - | tr \\t , | sed -e 's/(.*)//' -e 's/, /,/g' | \
		awk -F, '_[$1]{next}{t=$4}$4~/group/{t=$5}$5~/group/{t=$6}$2~/Virus|Metagenomic/{t=$2} {_[$1]=1;print$1,t}' OFS=, | \
		awk -F, 'NR==FNR{c[$1]=$3;next} c[$2]{print $1,c[$2],$2}' OFS=, datasets/taxonomy.list -)
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
	$(awk -F, 'NR==FNR{c[$1]=$3;next} c[$NF]{print $1,c[$NF],$NF}' OFS=, datasets/habitats.list habitats.csv)
EOF

cat <<- EOF > datasets/organism.txt
	DATASET_TEXT
	SEPARATOR COMMA
	DATASET_LABEL,Organism
	DATA
	$(tail -n+2 uniprot.tax | awk -F\\t '{NF--}1' OFS=\\t | cat uniprot.tax.override - | tr \\t , | sed -e 's/, /,/g' | awk -F, '{print$1,$NF}' OFS=, | \
		sed -e 's/(strain \(.*\))/\1/' -e 's/) (.*//g' -e 's/(.*//g' | awk -F, '!_[$1]{_[$1]=1; print$0,-1,"black","normal",1,0}' OFS=,)
EOF

cat <<- EOF > datasets/connections.txt
	DATASET_CONNECTION
	SEPARATOR COMMA
	DATASET_LABEL,Multiple NblA's
	DATA
	$(tail -n+2 uniprot.tax | awk -F\\t '{NF--}1' OFS=\\t | cat uniprot.tax.override - | tr \\t , | sed -e 's/, /,/g' | \
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