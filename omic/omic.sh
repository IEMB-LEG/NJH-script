

grep ">" .aa |sed  's/>//' > name
python omic-rename2.py blast2go.gff name out
python omic-gff-merge.py ../../intron/calkinsi-purged.gff3  b2 b2-merge.gff
