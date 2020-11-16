for i in $(seq 1 22) X Y M;
do echo $i;
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr${i}.fa.gz;
done
gunzip *.gz
for i in $(seq 1 22) X Y M;
do cat chr${i}.fa >> hg19.fasta;
done
#rm -fr chr*.fasta
