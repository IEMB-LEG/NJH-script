
genome=$1
gff3=$2
prefix=$3

python ~/script/intron/intron-ratio-stats.py $genome $gff3 $prefix-intron-length-ratio.txt
python ~/script/intron/intron-ratio-plot2.py $prefix-intron-length-ratio.txt intron-ratio-plot.pdf
python ~/script/intron/intron2.py $gff3 $prefix-intron.gff
python ~/script/intron/intron-extract2.py $genome $prefix-intron.gff $prefix-intron.fasta
python ~/script/intron/plot_contour_script3.py $prefix-intron.fasta $prefix-intron.stats
cut -f 2-4 $prefix-intron.stats |sed '1,8d'> $prefix-intron-length-coden.txt

python ~/script/intron/plot_cont5.py $prefix-intron-length-coden.txt $prefix-intron-length-coden-3dev.pdf --interval 10 --max_length 70 --y_max 100000
python ~/script/intron/plot_cont4.py $prefix-intron-length-coden.txt $prefix-intron-length-coden-no3dev.pdf --interval 10 --max_length 70 --y_max 100000
python intron-gc-content3.py $prefix-intron.stats $prefix-intron-gc.pdf
python calculate_codon_usage2.py $genome $gff3 $prefix-coden-usage.pdf
#python ~/script/intron/number-plot5.py alkinsi-purged.gff3 duboscqui-purged.gff3 nephridiatum-purged.gff3 putrinum-purged.gff3 woodruffi-purged.gff3 tetraurelia_mac_d4-2_annotation_v1.gff3 all2.plot.pdf
#python ~/script/intron/number-plot7.py --max_gene_length 5000 --max_exon_count 10 --max_intron_length 50 --max_exon_length 1400 calkinsi-purged.gff3 duboscqui-purged.gff3 nephridiatum-purged.gff3 putrinum-purged.gff3 woodruffi-purged.gff3 tetraurelia_mac_d4-2_annotation_v1.gff3 all3.plot.pdf
