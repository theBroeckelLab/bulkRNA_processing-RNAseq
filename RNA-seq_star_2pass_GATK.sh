#!/bin/bash

echo $1

start="$(date +%s)"

filename=$(basename $1)
extension=${filename##*.}
filename=${filename%.*}
sample=$"${filename}"
echo $sample
starSAMprefix=$"${filename}_STAR_pass2"
starSAM=$"${filename}_STAR_pass2Aligned.out.sam"
starRGaddedSortedBAM=$"${filename}_STAR_pass2_RG_sorted.bam"
starRGaddedSortedMarkdupBAM=$"${filename}_STAR_pass2_markdup.bam"
MarkdupMetrics=$"${filename}_STAR_pass2_markdup.metrics"
starMarkdupSplitBAM=$"${filename}_pass2_markdup_NCigarsplit.bam"
output_vcf=$"${filename}_GATK_haplotypecaller.vcf"
output_filtered_vcf=$"${filename}_GATK_haplotypecaller_filtered.vcf"
STARbowtieHTSeqGeneCount=$"${filename}_STAR_bowtie2_combined_htseq_ucsc_hg19_genes.counts"
STARbowtieHTSeqEnsembleCount=$"${filename}_STAR_bowtie2_combined_htseq_ensemble_75.counts"
STARbowtieHTSeqGencodeCount=$"${filename}_STAR_bowtie2_combined_htseq_gencode_v19.counts" 
STARbowtieHTSeqrRNACount=$"${filename}_STAR_bowtie2_combined_htseq_ucsc_hg19_rrna.counts"

STAR --genomeDir /data/IonTorrent/Proton/RNA-seq/mRNA_2014/myCells/test_run/processed_files/star_pass2/genomedir --readFilesIn $1 --runThreadN 10 --outFileNamePrefix $starSAMprefix --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outSAMstrandField intronMotif --chimSegmentMin 15 --alignIntronMax 50000 --seedSearchStartLmax 30

java -Xmx4g -jar /usr/local/picard/AddOrReplaceReadGroups.jar I=$starSAM O=$starRGaddedSortedBAM SO=coordinate RGID=$sample RGLB="RNA" RGPL="IONTORRENT" RGPU="Proton" RGSM=$sample VALIDATION_STRINGENCY=LENIENT

java -Xmx4g -jar /usr/local/picard/MarkDuplicates.jar I=$starRGaddedSortedBAM O=$starRGaddedSortedMarkdupBAM CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=$MarkdupMetrics

java -Xmx8g -jar /data/Non_system_tools/GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $starRGaddedSortedMarkdupBAM -o $starMarkdupSplitBAM -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

#### put in a command for the base recalibration

java -Xmx8g -jar /data/Non_system_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 4 -R /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $starMarkdupSplitBAM -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $output_vcf

java -Xmx8g -jar /data/Non_system_tools/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -V $output_vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $output_filtered_vcf

htseq-count -q -i gene_name $starSAM /ref/bedtools_refs/gencode.v19.annotation.gtf > $STARbowtieHTSeqGencodeCount
htseq-count -q -i gene_id $starSAM /ref/bedtools_refs/rrna_hg19.gtf > $STARbowtieHTSeqrRNACount

finish="$(date +%s)"
echo "Analysis complete! Total time in seconds : $((finish-start))" | mail -s "RNA-Seq analysis using STAR 2-pass complete for file $1" paggarwal@mcw.edu
