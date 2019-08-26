#!/bin/bash
DIR=/home/chenzh/My_project/NG_smallRNA
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BDOC=$DIR/big_doc
BIN=$DIR/bin

cd $DIR
## bash smallrna_work_driver.bash D10_14 

SA=$1
WD=$TMPD/raw_data/${SA}
cutada=$BDOC/cutadapt_3prime.fa ### From paper
RGI=$BDOC/RefSeq/star_index
REF=$BDOC/RefSeq/genome.fa

cd $WD

umi_tools extract --bc-pattern=NNNNNNNN -I ${SA}.fastq -S ${SA}.extract.fastq -L ${SA}.extract.log

cutadapt -a file:${cutada}   -e 0.1 -O 1 -m 18 -u 2 -o ${SA}.extract.cutadp.fastq ${SA}.extract.fastq  >${SA}.extract.cutadp.log


STAR --runThreadN 2 --genomeDir ${RGI} --readFilesIn ${SA}.extract.cutadp.fastq --outSAMstrandField intronMotif --outFilterMultimapNmax 50   --outFilterScoreMinOverLread 0   --outFilterMatchNmin 18   --outFilterMatchNminOverLread 0   --outFilterMismatchNoverLmax 0.04   --alignIntronMax 1 --readFilesCommand - --outSAMtype BAM SortedByCoordinate  --outFileNamePrefix ${SA}.


python $SRC/remove_softclipped_reads_single.py -f 0 -t 3 -m 40 -i ${SA}.Aligned.sortedByCoord.out.bam -w ${SA}.remove.sc.bam -o ${SA}.remove.sc.max40.bam > ${SA}.rsclp.log

samtools index ${SA}.remove.sc.max40.bam
samtools index ${SA}.remove.sc.bam

#### remove duplication
umi_tools dedup --method adjacency --output-stats dedup.log1  --read-length -I ${SA}.remove.sc.max40.bam -S ${SA}.remove.sc.max40.dedup.bam > /dev/null
umi_tools dedup --method adjacency --output-stats dedup.log2  --read-length -I ${SA}.remove.sc.bam -S ${SA}.remove.sc.dedup.bam > /dev/null

### remove false positive adapter
python $SRC/remove_reads_with_genomic_AGG_single.py -g $REF -x 41 -i ${SA}.remove.sc.max40.dedup.bam -o ${SA}.remove.sc.max40.dedup.AGG.bam
####

#### Get counts
python $SRC/count_smallrnas.modified.py -g $BDOC/small.anno.bed -i ${SA}.remove.sc.max40.dedup.AGG.bam -o ${SA}
python $SRC/count_smallrnas.modified.py -g $BDOC/prec.anno.bed -i ${SA}.remove.sc.dedup.bam -o ${SA}.pre



echo -e "#samples\t"$SA  > ${SA}.mature.count.out
grep -v "^#" ${SA}_Count.txt |sort -k1,1  |groupBy -g 1 -c 3 -o sum|perl $SRC/become.format.pl $BDOC/mature.out.format  >> ${SA}.mature.count.out

echo -e "#samples\t"$SA  > ${SA}.prec.count.out
grep -v "^#" ${SA}.pre_Count.txt |sort -k1,1  |groupBy -g 1 -c 3 -o sum|perl $SRC/become.format.pl $BDOC/prec.out.format  >> ${SA}.prec.count.out

