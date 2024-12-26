# mkdir telomere
# ls *.fasta > list2
# cd telomere
# while read a ; do bash ~/script/telomere.sh $a ; done < ../list2
# contig should be longer than 300, extract first and last 300nt from each contig, search "TTGGGGTTGGGG" & "CCCCAACCCCAA" and sum.
#!/bin/bash
i=$1
mkdir telomere
cd telomere
awk '{if($0 ~ />/) {print $0} else {printf("%s", $0)}}' ../$i | sed 's/>/\n>/g' | sed '$ s/$/\n/'| sed 1d > n1$i
awk -vn=300 '{print substr($0,length($0)-n+1)}' n1$i >n1-last300$i
awk '{print substr($0,1,300)}' n1$i > n1-first300$i
paste n1-first300$i n1-last300$i > n1-first-last300$i
sed -i 's/>.*>/>/' n1-first-last300$i
sed -i 's/	//' n1-first-last300$i
grep "TTTGGGTTTGGG" -B1 n1-first-last300$i > $i-tttggg.txt
grep "CCCAAACCCAAA" -B1 n1-first-last300$i > $i-cccaaa.txt
grep ">" $i-tttggg.txt > $i-tttggg-contig.txt
grep ">" $i-cccaaa.txt > $i-cccaaa-contig.txt
wc -l $i-cccaaa-contig.txt
wc -l $i-tttggg-contig.txt
grep -wf $i-tttggg-contig.txt $i-cccaaa-contig.txt | wc -l
grep "TTGGGGTTGGGG" -B1 n1-first-last300$i > $i-ttgggg.txt
grep "CCCCAACCCCAA" -B1 n1-first-last300$i > $i-ccccaa.txt
grep ">" $i-ttgggg.txt > $i-ttgggg-contig.txt
grep ">" $i-ccccaa.txt > $i-ccccaa-contig.txt
wc -l $i-ccccaa-contig.txt
wc -l $i-ttgggg-contig.txt
grep -wf $i-ttgggg-contig.txt $i-ccccaa-contig.txt | wc -l
cat $i-tttggg-contig.txt $i-ttgggg-contig.txt |sort |uniq > $i-tg-contig.txt
cat $i-cccaaa-contig.txt $i-ccccaa-contig.txt |sort |uniq > $i-ca-contig.txt
echo "tttggg & ttgggg:"
wc -l $i-tg-contig.txt
echo "cccaaa & ccccaa:"
wc -l $i-ca-contig.txt
echo "complete:"
grep -wf $i-tg-contig.txt $i-ca-contig.txt | wc -l
#cat $i-ttgggg-contig.txt $i-ccccaa-contig.txt |sort |uniq -d|wc -l

