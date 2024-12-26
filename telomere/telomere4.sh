# while read a ; do bash ~/script/telomere4.sh $a ; done < list3
# use longer sequence for cd hit, exam telomere2 clustered contigs.
i=$1
cd $i-cdhit
grep "at" -B1 $i-ttgggg-cdhit.fasta.clstr | cut -f 2 -d " " > $i-tg.contig
sed -i 's/>//' $i-tg.contig
sed -i 's/\..*//' $i-tg.contig
grep "at" -B1 $i-ccccaa-cdhit.fasta.clstr | cut -f 2 -d " " > $i-ca.contig
sed -i 's/>//' $i-ca.contig
sed -i 's/\..*//' $i-ca.contig
seqkit grep -j 6  --pattern-file $i-tg.contig ~/data/yue-hifi/genome/catgenome/$i-cat.fasta > $i-tg-cluster.fasta
seqkit grep -j 6  --pattern-file $i-ca.contig ~/data/yue-hifi/genome/catgenome/$i-cat.fasta > $i-ca-cluster.fasta
awk '{if($0 ~ />/) {print $0} else {printf("%s", $0)}}' $i-ca-cluster.fasta | sed 's/>/\n>/g' | sed '$ s/$/\n/'| sed 1d > n1-$i-ca-cluster.fasta
awk '{if($0 ~ />/) {print $0} else {printf("%s", $0)}}' $i-tg-cluster.fasta | sed 's/>/\n>/g' | sed '$ s/$/\n/'| sed 1d > n1-$i-tg-cluster.fasta
awk -vn=100000 '{print substr($0,length($0)-n+1)}' n1-$i-tg-cluster.fasta > n1-$i-tg-cluster-last.fasta
awk '{print substr($0,1,100000)}' n1-$i-ca-cluster.fasta > n1-$i-ca-cluster-first.fasta

cdhit -i n1-$i-tg-cluster-last.fasta -o n1-$i-tg-cluster-last-exam.fasta -c 0.85 -d 0 -aL 0 -AL 100000 -aS 0.9 -AS 10000 -T 8 -M 24000
cdhit -i n1-$i-ca-cluster-first.fasta -o n1-$i-ca-cluster-first-exam.fasta -c 0.85 -d 0 -aL 0 -AL 100000 -aS 0.9 -AS 10000 -T 8 -M 24000

grep "*"  n1-$i-ca-cluster-first-exam.fasta.clstr | cut -f 2 -d " " > $i-ca-telo.contig
grep "*" $i-ttgggg-cdhit.fasta.clstr | cut -f 2 -d " " >> $i-ca-telo.contig
sort $i-ca-telo.contig | uniq > $i-ca-telomere.contig
sed -i 's/>//' $i-ca-telomere.contig
sed -i 's/\..*//' $i-ca-telomere.contig
grep "*" n1-$i-tg-cluster-last-exam.fasta.clstr | cut -f 2 -d " " > $i-tg-telo.contig
grep "*" $i-ccccaa-cdhit.fasta.clstr | cut -f 2 -d " " >> $i-tg-telo.contig
sort $i-tg-telo.contig | uniq > $i-tg-telomere.contig
sed -i 's/>//' $i-tg-telomere.contig
sed -i 's/\..*//' $i-tg-telomere.contig

cd ..
