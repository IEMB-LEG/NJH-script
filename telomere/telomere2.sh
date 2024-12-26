# while read a ; do bash ~/script/telomere2.sh $a ; done < ../list3
i=$1
mkdir subreads-$i
cd subreads-$i

sed -i 's/>//' ../$i-cat.fasta-ccccaa-contig.txt
sed -i 's/>//' ../$i-cat.fasta-ttgggg-contig.txt

seqkit grep -j 6  --pattern-file ../$i-cat.fasta-ccccaa-contig.txt ../../$i-cat.fasta > $i.ccccaa.sub.fasta
seqkit sort $i.ccccaa.sub.fasta > $i.ccccaa.subreads.fasta
seqkit stats -j 6 $i.ccccaa.subreads.fasta

awk '{if($0 ~ />/) {print $0} else {printf("%s", $0)}}' $i.ccccaa.subreads.fasta | sed 's/>/\n>/g' | sed '$ s/$/\n/'| sed 1d > n1$i.ccccaa.subreads.fasta
awk '{print substr($0,1,10000)}' n1$i.ccccaa.subreads.fasta > n1-first10000-$i.ccccaa.subreads.fasta

seqkit grep -j 6  --pattern-file ../$i-cat.fasta-ttgggg-contig.txt ../../$i-cat.fasta > $i.ttgggg.sub.fasta
seqkit sort $i.ttgggg.sub.fasta > $i.ttgggg.subreads.fasta
seqkit stats -j 6 $i.ttgggg.subreads.fasta
awk '{if($0 ~ />/) {print $0} else {printf("%s", $0)}}' $i.ttgggg.subreads.fasta | sed 's/>/\n>/g' | sed '$ s/$/\n/'| sed 1d > n1$i.ttgggg.subreads.fasta
awk -vn=10000 '{print substr($0,length($0)-n+1)}' n1$i.ttgggg.subreads.fasta >n1-last10000-$i.ttgggg.subreads.fasta

cd ..
