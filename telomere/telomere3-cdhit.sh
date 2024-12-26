# while read a ; do bash ~/script/telomere3-cdhit.sh $a ; done < ../list3
i=$1
mkdir $i-cdhit
cd $i-cdhit
cdhit -i ../subreads-$i/n1-first10000-$i.ccccaa.subreads.fasta -o $i-ccccaa-cdhit.fasta -c 0.85 -d 0 -aL 0.9 -AL 1000 -aS 0.9 -AS 1000 -T 8 -M 24000
cdhit -i ../subreads-$i/n1-last10000-$i.ttgggg.subreads.fasta -o $i-ttgggg-cdhit.fasta -c 0.85 -d 0 -aL 0.9 -AL 1000 -aS 0.9 -AS 1000 -T 8 -M 24000
cd ..
#-i：输入文件，fasta格式。
#-o：输出文件前缀，输出文件有两个，分别为fasta格式序列文件和以.clstr结尾的聚类信息文件。
#-c：较短序列比对到长序列的bp与自身bp数的比值超过该数值则聚类为一组，默认为0.9。
#-d：聚类信息文件中各个聚类组中序列名的长度，设为0则将取完整序列名。
#-aL：控制代表序列比对严格程度的参数，默认为0，若设为0.8则表示比对区间要占到代表（长）序列的80%。
#-AL：控制代表序列比对严格程度的参数，默认为99999999，若设为40则表示代表序列的非比对区间要短于40bp。
#-aS：控制短序列比对严格程度的参数，默认为0，若设为0.8则表示比对区间要占到短序列的80%。
#-AS：控制短序列比对严格程度的参数，默认为99999999，若设为40则表示短序列的非比对区间要短于40bp。
####c is important, if length is similar, aL high; length is different, aL low, aS high, for more accuracy
