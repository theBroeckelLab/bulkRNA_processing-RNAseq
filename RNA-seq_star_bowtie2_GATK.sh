#!/bin/bash -e

echo $1

## Starting time

start="$(date +%s)"

## output usage if only the script is provided and no input file is given

if test "$1" == ""; then
	echo $'\a'Usage: /data/Non_system_tools/scripts/rna_seq_scripts/RNA-seq_star_bowtie2_GATK.sh fastq_file
	exit
fi

## declare all the variables

filename=$(basename $1)
extension=${filename##*.}
filename=${filename%.*}
sample=$"${filename}"
starSAMprefix=$"${filename}_STAR_hg19"
starSAM=$"${filename}_STAR_hg19Aligned.out.sam"
##starBAM=$"${filename}_STAR_hg19.bam"
##starUnmappedSAM=$"${filename}_STAR_unmapped_hg19.sam"
starUnmappedFastq=$"${filename}_STAR_hg19Unmapped.out.mate1"
bowtieSAM=$"${filename}_STAR_unmapped_bowtie2_local_hg19.sam"
##bowtieBAM=$"${filename}_STAR_unmapped_bowtie2_local_hg19.bam"
bowtieMetrics=$"${filename}_STAR_unmapped_bowtie2_local_hg19.metrics"
STARbowtieSAM=$"${filename}_STAR_bowtie2_local_combined_hg19.sam"
RGaddedSortedBAM=$"${filename}_STAR_bowtie2_combined_readgroup.bam"
RGaddedSortedMarkdupBAM=$"${filename}_STAR_bowtie2_combined_markdup.bam"
##RGaddedSortedMarkdupMappedBAM=$"${filename}_STAR_bowtie2_combined_markdup_mapped.bam"
MarkdupMetrics=$"${filename}_markdup.metrics"
MarkdupSplitBAM=$"${filename}_markdup_cigarsplit.bam"
output_vcf=$"${filename}_GATK_haplotypecaller.vcf"
output_filtered_vcf=$"${filename}_GATK_haplotypecaller_filtered.vcf"
STARbowtieHTSeqGeneCount=$"${filename}_STAR_bowtie2_combined_htseq_ucsc_hg19_genes.counts"
##STARbowtieHTSeqEnsembleCount=$"${filename}_STAR_bowtie2_combined_htseq_ensemble_75.counts"
STARbowtieHTSeqGencodeCountsGT200=$"${filename}_STAR_bowtie2_combined_htseq_gencode_v19_gt200.counts"
STARbowtieHTSeqGencodeCountsGT200_reverse=$"${filename}_STAR_bowtie2_combined_htseq_gencode_v19_gt200_reverse.counts"
STARbowtieHTSeqGencodeCountsGT200Q20=$"${filename}_STAR_bowtie2_combined_htseq_gencode_v19_gt200_Q20.counts"
STARbowtieHTSeqGencodeCount=$"${filename}_STAR_bowtie2_combined_htseq_gencode_v19.counts"
STARbowtieHTSeqrRNACount=$"${filename}_STAR_bowtie2_combined_htseq_ucsc_hg19_rrna.counts"
IntermediateDirectory=$"${filename}_temp_dir"

if [ -d $IntermediateDirectory ]; then
   echo $'\a'Error: temp directory already exists. Check it!
   exit
fi

mkdir $IntermediateDirectory

#### Start Analysis

## Run FastQC

#start_fastqc="$(date +%s)"

#/data/Non_system_tools/FastQC/fastqc $1

#finish_fastqc="$(date +%s)"
#echo "Time taken by FastQC, in seconds : $((finish_fastqc-start_fastqc))" | mail -s "FastQC for $sample complete" paggarwal@mcw.edu

## STAR alignment

start_star="$(date +%s)"

STAR --genomeDir /ref/STAR_ref/STARIndex_ucsc_hg19_gtf/ --readFilesIn $1 --runThreadN 10 --outFileNamePrefix $starSAMprefix --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --chimSegmentMin 15 --alignIntronMax 50000 --seedSearchStartLmax 30

#STAR --genomeDir /ref/STAR_ref/STARIndex_ucsc_hg19_gtf/ --readFilesIn $1 --runThreadN 10 --outFileNamePrefix $starSAMprefix --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.04 --outReadsUnmapped Fastx --chimSegmentMin 15 --alignIntronMax 50000 --seedSearchStartLmax 30

finish_star="$(date +%s)"
echo "Time taken by STAR, in seconds : $((finish_star-start_star))" | mail -s "STAR alignment for $sample complete" paggarwal@mcw.edu

## Bowtie2 local alignment

start_bowtie2="$(date +%s)"

/usr/local/bowtie2-2.2.2/bowtie2 --very-sensitive-local -p 10 -x /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -U $starUnmappedFastq -S $bowtieSAM --no-unal --met-file $bowtieMetrics

finish_bowtie2="$(date +%s)"
echo "Time taken to run Bowtie2, in seconds : $((finish_bowtie2-start_bowtie2))" | mail -s "Bowtie2 alignment for $sample complete" paggarwal@mcw.edu

## Merge STAR and Bowtie2 alignments

start_merge="$(date +%s)"

java -Xmx8g -jar /usr/local/picard/MergeSamFiles.jar INPUT=$starSAM INPUT=$bowtieSAM OUTPUT=$STARbowtieSAM VALIDATION_STRINGENCY=LENIENT

finish_merge="$(date +%s)"
echo "Time taken in merging the SAM files, in seconds : $((finish_merge-start_merge))" | mail -s "Merging alignments for $sample complete" paggarwal@mcw.edu

## Create read group markdup BAM files

start_markdup="$(date +%s)"

java -Xmx8g -jar /usr/local/picard/AddOrReplaceReadGroups.jar I=$STARbowtieSAM O=$RGaddedSortedBAM SO=coordinate RGID=$sample RGLB="RNA" RGPL="IONTORRENT" RGPU="Proton" RGSM=$sample VALIDATION_STRINGENCY=LENIENT

java -Xmx4g -jar /usr/local/picard/MarkDuplicates.jar I=$RGaddedSortedBAM O=$RGaddedSortedMarkdupBAM CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT M=$MarkdupMetrics

finish_markdup="$(date +%s)"
echo "Time taken to create ReadGroup/Markdup BAM files, in seconds : $((finish_markdup-start_markdup))" | mail -s "ReadGroup/Markdup BAM created for $sample" paggarwal@mcw.edu

## HTSeq quantification

start_htseq="$(date +%s)"
##htseq-count -q -i gene_id $STARbowtieSAM /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf > $STARbowtieHTSeqGeneCount
##htseq-count -q -i gene_name $STARbowtieSAM /ref/bedtools_refs/Homo_sapiens.GRCh37.75.edited.gtf > $STARbowtieHTSeqEnsembleCount
htseq-count -q -i gene_name $STARbowtieSAM /ref/bedtools_refs/gencode_v19_annotation_wo_genes_smaller_than_200bp.gtf > $STARbowtieHTSeqGencodeCountsGT200
#htseq-count -q -s reverse -i gene_name $STARbowtieSAM /ref/bedtools_refs/gencode_v19_annotation_wo_genes_smaller_than_200bp.gtf > $STARbowtieHTSeqGencodeCountsGT200_reverse
htseq-count -q -i gene_name $STARbowtieSAM /ref/bedtools_refs/gencode.v19.annotation.gtf > $STARbowtieHTSeqGencodeCount
##htseq-count -q -i gene_name $STARbowtieSAM /ref/bedtools_refs/gencode_v19_annotation_wo_genes_smaller_than_200bp.gtf > $STARbowtieHTSeqGencodeCountsGT200
htseq-count -q -i gene_name -a 20 $STARbowtieSAM /ref/bedtools_refs/gencode_v19_annotation_wo_genes_smaller_than_200bp.gtf > $STARbowtieHTSeqGencodeCountsGT200Q20
htseq-count -q -i gene_id $STARbowtieSAM /ref/bedtools_refs/rrna_hg19.gtf > $STARbowtieHTSeqrRNACount

finish_htseq="$(date +%s)"
echo "Time taken for HTSeq, in seconds : $((finish_htseq-start_htseq))" | mail -s "HTSeq quanitification for $sample complete" paggarwal@mcw.edu

## Prepare aligned file for variant calling

start_prep="$(date +%s)"

java -Xmx8g -jar /data/Non_system_tools/GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $RGaddedSortedMarkdupBAM -o $MarkdupSplitBAM -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

finish_prep="$(date +%s)"
echo "Time taken to prepare files, in seconds : $((finish_prep-start_prep))" | mail -s "File preparation for variant calling for $sample complete" paggarwal@mcw.edu

#### put in a command for the base recalibration

## GATK HaplotypeCaller

start_hap="$(date +%s)"

java -Xmx8g -jar /data/Non_system_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 4 -R /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -I $MarkdupSplitBAM -recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $output_vcf

java -Xmx8g -jar /data/Non_system_tools/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R /ref/UCSC_hg19_2014/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -V $output_vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $output_filtered_vcf

finish_hap="$(date +%s)"
echo "Time taken for GATK calling and filtering, in seconds : $((finish_hap-start_hap))" | mail -s "GATK based variant calling for $sample complete" paggarwal@mcw.edu

## Run FastQC

start_fastqc="$(date +%s)"

/data/Non_system_tools/FastQC/fastqc $1

finish_fastqc="$(date +%s)"
echo "Time taken by FastQC, in seconds : $((finish_fastqc-start_fastqc))" | mail -s "FastQC for $sample complete" paggarwal@mcw.edu

### Analysis complete

finish="$(date +%s)"
echo "Analysis complete! Total time, in seconds : $((finish-start))" | mail -s "RNA-Seq analysis using Bowtie2 complete for file $sample" paggarwal@mcw.edu

## Move all the temp files to the temp directory

mv $starSAM $IntermediateDirectory
mv $bowtieSAM $IntermediateDirectory
mv $MarkdupSplitBAM $IntermediateDirectory
mv $STARbowtieSAM $IntermediateDirectory
mv $starUnmappedFastq $IntermediateDirectory
mv $output_vcf $IntermediateDirectory
