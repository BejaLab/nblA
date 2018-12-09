
json=$(realpath "$1")
svg=$2
awk -v json="$json" 'BEGIN{_="jsonFile=\"" json "\";\n"}NR==FNR{_=_$0"\n";next}1;/<svg/{print "<script type=\"text/javascript\"><![CDATA["; print _; print"]]></script>"}' skylign.js "$svg"
