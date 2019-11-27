#!/usr/bin/env bash
#SBATCH -A SNIC2019-4-10
#SBATCH -p hebbe
#SBATCH -J aligner
#SBATCH -n 2
#SBATCH -t 0-08:00:00

module load Anaconda3/5.3.0;
source activate /c3se/NOBACKUP/users/mariaji/envs/bwa;
module load iccifort/2019.1.144-GCC-8.2.0-2.31.1;
module load SAMtools/1.9;
module load GATK/3.7-Java-1.8.0_112;
module load picard/2.18.25-Java-1.8;

# bwa index -a bwtsw ${REF};
# samtools faidx ${REF}

# java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
# REFERENCE=${REF} \

# mkdir rawbams;
# mkdir indelvcf;
# mkdir realnbam;
# mkdir rawvcf;
# mkdir recal1;
# mkdir recal2;
# mkdir variants;

#######################################################################
# 1. map nodup reads to loci reference abd prepare bam
#######################################################################

# modify depending on data
# running from outside folder with index folders
# reads are trimmed and PCR duplicates masked with Picard
SAMPLE=$(ls Index${SLURM_ARRAY_TASK_ID}/*_nodup_clean-READ1.fastq | sed -r 's/_nodup_clean-READ1.fastq//g');
GROUP=$(echo ${SAMPLE} | sed -r 's/Index[0-9]+\///g')
REF=$(path/to/ref.fasta)
# echo ${SAMPLE} >> text.txt
# echo ${GROUP} >> text.txt

bwa mem -t 4 ${REF} \
    ${SAMPLE}_nodup_clean-READ1.fastq \
    ${SAMPLE}_nodup_clean-READ2.fastq \
    -o rawbams/${SAMPLE}.sam $2> rawbams/${SAMPLE}_step1.out;

java -jar $EBROOTPICARD/picard.jar SamFormatConverter \
    I=rawbams/${SAMPLE}.sam \
    O=rawbams/${SAMPLE}.bam;

java -jar $EBROOTPICARD/picard.jar SortSam \
      I=rawbams/${SAMPLE}.bam \
      O=rawbams/${SAMPLE}.sorted.bam \
      SORT_ORDER=coordinate

java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
      I=rawbams/${SAMPLE}.sorted.bam;

# when automatising this make sure that the ${SAMPLE} variable excludes folder paths
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
       I=rawbams/${SAMPLE}.sorted.bam \
       O=rawbams/${SAMPLE}.sorted.wgroups.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=${GROUP} $2> rawbams/${SAMPLE}_groupfix.out;

java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
      I=rawbams/${SAMPLE}.sorted.wgroups.bam;

rm rawbams/${SAMPLE}.sorted.bam;
rm rawbams/${SAMPLE}.sorted.bai;
rm rawbams/${SAMPLE}.sam;
rm rawbams/${SAMPLE}.bam;

echo '**** Index built and bam-sam-bam switch performed****';

#######################################################################
# 2. call first VCF
#######################################################################

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T UnifiedGenotyper \
   -R ${REF} \
   -I rawbams/${SAMPLE}.sorted.wgroups.bam \
   -o rawvcf/${SAMPLE}.raw.vcf $2> rawvcf/${SAMPLE}_step2.out;

echo '****raw variants called****';

######################################################################
# 3. re-align indels
######################################################################

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R ${REF} \
   -V rawvcf/${SAMPLE}.raw.vcf \
   -o indelvcf/${SAMPLE}_indelraw.vcf \
   -selectType INDEL \
   --minIndelSize 2 \
   --maxIndelSize 50 $2> indelvcf/${SAMPLE}_step3_callindel.out;

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R ${REF} \
   -I rawbams/${SAMPLE}.sorted.wgroups.bam \
   --known indelvcf/${SAMPLE}_indelraw.vcf \
   -o indelvcf/${SAMPLE}.indel.intervals $2> indelvcf/${SAMPLE}_step3_indelinterv.out;

# targetIntervals output from RealignerTargetCreator
java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R ${REF} \
   -I rawbams/${SAMPLE}.sorted.wgroups.bam \
   -known indelvcf/${SAMPLE}_indelraw.vcf \
   -targetIntervals indelvcf/${SAMPLE}.indel.intervals \
   -o realnbam/${SAMPLE}.realn.bam $2> indelvcf/${SAMPLE}_step3_realn.out;
# need to deal with "Not attempting realignment in interval ... because there are too many reads."

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T UnifiedGenotyper \
   -R ${REF} \
   -I realnbam/${SAMPLE}.realn.bam \
   -o realnbam/${SAMPLE}.realn.vcf $2> realnbam/${SAMPLE}_step3_ralnvcf.out;

rm rawbams/${SAMPLE}.sorted.wgroups.bam;
rm rawbams/${SAMPLE}.sorted.wgroups.bai;
rm indelvcf/${SAMPLE}_indelraw.vcf;
rm indelvcf/${SAMPLE}.indel.intervals;

#######################################################################
# 4. recalibration round 1 - ideally run two minimum
#######################################################################

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ${REF} \
    -I realnbam/${SAMPLE}.realn.bam \
    -knownSites realnbam/${SAMPLE}.realn.vcf \
    -o recal1/${SAMPLE}_recal_data.1.table $2> recal1/${SAMPLE}_step4_recal1.out;

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ${REF} \
    -I realnbam/${SAMPLE}.realn.bam \
    -knownSites realnbam/${SAMPLE}.realn.vcf \
    -BQSR recal1/${SAMPLE}_recal_data.1.table \
    -o recal1/${SAMPLE}_post_recal_data.1.table $2> recal1/${SAMPLE}_step4_recal1post.out; 
    
#######################################################################
# 5. filter with recalibration and call gVCF
#######################################################################

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T PrintReads \
   -R ${REF} \
   -I realnbam/${SAMPLE}.realn.bam \
   -BQSR recal1/${SAMPLE}_post_recal_data.1.table \
   -o recal1/${SAMPLE}.recab.bam $2> recal1/${SAMPLE}_step4_recab.out;

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -T UnifiedGenotyper \
   -R ${REF} \
   -I recal1/${SAMPLE}.recab.bam \
   -o variants/${SAMPLE}.nodup.recab.realn.vcf \
   --output_mode EMIT_ALL_SITES $2> variants/${SAMPLE}_varcal.out;

rm realnbam/${SAMPLE}.realn.bam

# rm realnbam/${SAMPLE}_realn.bam;
touch "${SAMPLE}"_done.out;
