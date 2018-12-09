awk -F\\t 'NR==FNR{_[">"$1]=1}/>/{p=0}/>/&&_[$1]{p=1}p' alignment-choice.list ../nbla.mafft.fasta > alignment-choice-selected.fasta
cols=$(trimal -sgc -in alignment-choice-selected.fasta | awk '$2==100{print$1}' ORS=,)
trimal -select { $cols } -in ../nbla.mafft.fasta | tr ' ' \\t | awk -F\\t 'NR==FNR{_[">"$1]=$2;next}/>/{p=0}/>/&&_[$1]{$0=">"_[$1];p=1}!p{l[NR]=$0}p;END{for(i=1;i<=NR;i++)if(l[i])print l[i]}' alignment-choice.list - > alignment-choice.fasta
