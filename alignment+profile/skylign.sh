awk 'BEGIN{_="jsonFile=\"/home/har-wradim/Documents/omer/nblA/alignment+profile/skylign.json\";\n"}NR==FNR{_=_$0"\n";next}1;/<svg/{print "<script type=\"text/javascript\"><![CDATA["; print _; print"]]></script>"}' skylign.js skylign.svg > skylign-clustalX.svg
