#!/bin/bash

echo $1

## Starting time

start="$(date +%s)"


filename=$(basename $1)
extension=${filename##*.}
filename=${filename%.*}
starSAMprefix=$"${filename}_STAR_hg19"
starSAM=$"${filename}_STAR_hg19Aligned.out.sam"
starBAM=$"${filename}_STAR_hg19.bam"
starUnmappedSAM=$"${filename}_STAR_unmapped_hg19.sam"
starUnmappedFastq=$"${filename}_STAR_unmapped_hg19.fastq"
bowtieSAM=$"${filename}_STAR_unmapped_bowtie2_local_hg19.sam"
bowtieBAM=$"${filename}_STAR_unmapped_bowtie2_local_hg19.bam"
bowtieMetrics=$"${filename}_STAR_unmapped_bowtie2_local_hg19.metrics"
STARbowtieSAM=$"${filename}_STAR_bowtie2_local_combined_hg19.sam"
STARbowtieBAM=$"${filename}_STAR_bowtie2_local_combined_hg19.bam"
STARbowtieHTSeqGeneCount=$"${filename}_STAR_bowtie2_combined_htseq_ucsc_hg19_genes.counts"
STARbowtieHTSeqEnsembleCount=$"${filename}_STAR_bowtie2_combined_htseq_ensemble_75.counts"
STARbowtieHTSeqGencodeCount=$"${filename}_STAR_bowtie2_combined_htseq_gencode_v19.counts"
STARbowtieHTSeqrRNACount=$"${filename}_STAR_bowtie2_combined_htseq_ucsc_hg19_rrna.counts"

STAR --genomeDir /ref/STAR_ref/STARIndex_ucsc_hg19_gtf/ --readFilesIn $1 --runThreadN 10 --outFileNamePrefix $starSAMprefix --outFilterMultimapNmax 1 --outFilterMismatchNmax 5 --outSAMstrandField intronMotif --outSAMunmapped Within --chimSegmentMin 15 --alignIntronMax 50000 --seedSearchStartLmax 30

samtools view -bS $starSAM > $starBAM

samtools view -f 4 $starBAM > $starUnmappedSAM

java -Xmx4g -jar /usr/local/picard/SamToFastq.jar INPUT=$starUnmappedSAM FASTQ=$starUnmappedFastq VALIDATION_STRINGENCY=LENIENT

/usr/local/bowtie2-2.2.2/bowtie2 --very-sensitive-local -p 4 -x /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -U $starUnmappedFastq -S $bowtieSAM --met-file $bowtieMetrics

java -Xmx4g -jar /usr/local/picard/MergeSamFiles.jar INPUT=$starSAM INPUT=$bowtieSAM OUTPUT=$STARbowtieSAM VALIDATION_STRINGENCY=LENIENT

samtools view -bS $STARbowtieSAM > $STARbowtieBAM

htseq-count -q -i gene_id $STARbowtieSAM /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > $STARbowtieHTSeqGeneCount
htseq-count -q -i gene_name $STARbowtieSAM /ref/bedtools_refs/Homo_sapiens.GRCh37.75.edited.gtf > $STARbowtieHTSeqEnsembleCount
htseq-count -q -i gene_name $STARbowtieSAM /ref/bedtools_refs/gencode.v19.annotation.gtf > $STARbowtieHTSeqGencodeCount
htseq-count -q -i gene_id $STARbowtieSAM /ref/bedtools_refs/rrna_hg19.gtf > $STARbowtieHTSeqrRNACount

finish="$(date +%s)"
echo "Analysis complete! Total time in seconds : $((finish-start))" | mail -s "RNA-Seq analysis complete for file $1" paggarwal@mcw.edu
