# Extract telomere contigs from purged genome. add telomere contigs from unpurged genome to purged genome
# while read a ; do bash ~/script/telomere5.sh $a
#cp /root/data/yue-hifi/genome/genome-final_hicanu_purged/gc40filter-1/*sort.fasta /root/data/yue-hifi/genome/catgenome/
i=$1
mkdir $i-exact
cd $i-exact
grep ">" ../../$i\gc40.sort.fasta > $i-purge.contig
cat ../../../genome-final_hicanu_purged/telomere/$i\gc40.sort.fasta-ccccaa-contig.txt ../../../genome-final_hicanu_purged/telomere/$i\gc40.sort.fasta-ttgggg-contig.txt | sort | uniq > $i-purge-telomere.contig
grep -vf $i-purge-telomere.contig $i-purge.contig > $i-purge.extract.contig
sed -i 's/>//' $i-purge.extract.contig
seqkit grep -j 8 --pattern-file $i-purge.extract.contig ../../$i\gc40.sort.fasta > $i-purge.extr.fasta
seqkit sort $i-purge.extr.fasta > $i-purge.extract.fasta
cd ../$i-cdhit
cat $i-ca-telomere.contig $i-tg-telomere.contig | sort | uniq  > $i-cat-telomere.contig
seqkit grep -j 8 --pattern-file $i-cat-telomere.contig ../../$i-cat.fasta > $i-cat-extr.fasta
seqkit sort $i-cat-extr.fasta > $i-cat.extract.fasta
cat ../$i-exact/$i-purge.extract.fasta $i-cat.extract.fasta > $i-gc40-telo.fasta
seqkit sort $i-gc40-telo.fasta > ../../$i-gc40-telomere.fasta


cd ..



