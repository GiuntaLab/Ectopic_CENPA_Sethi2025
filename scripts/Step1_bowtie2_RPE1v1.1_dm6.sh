#!/bin/sh
#Script for alignment and read count steps

#Indexed file path and mapping parameters
Idx='/data/index_reference/RPE1v1.1_dm6.fasta'
Set='--no-mixed --no-discordant --no-unal'

#Creating output directories if they don't already exist
if [ ! -d  Bam ];then mkdir Bam; fi
if [ ! -d  Log ];then mkdir Log; fi
if [ ! -d  ReadCount ];then mkdir ReadCount; fi

#Loop over all files matching the pattern
for k in `ls *_R1*.fastq.gz`

do

#Defining the variable i which serves as the base sample name
i=${k%_R1*.fastq.gz}

#Mapping
bowtie2 -p 15 -x $Idx $Set \
-1 ${i}_R1*.fastq.gz \
-2 ${i}_R2*.fastq.gz \
2>> Log/Log_Bowtie2_${i}.txt \
| samtools view -bS - | samtools sort - -@ 8 -O BAM -o ${i}_RPE1v1.1_dm6_sort.bam

#Indexing
samtools index -@ 8 ${i}_RPE1v1.1_dm6_sort.bam

#Extraction of RPE1 mapped reads per haplotype
samtools view -b ${i}_RPE1v1.1_dm6_sort.bam chr{1..22}_hap1 chrX_hap1 > ${i}_RPE1v1.1_hap1_sort.bam
samtools view -b ${i}_RPE1v1.1_dm6_sort.bam chr{1..22}_hap2 chrX_hap2 > ${i}_RPE1v1.1_hap2_sort.bam

#Extraction of dm6 mapped reads
samtools view -b ${i}_RPE1v1.1_dm6_sort.bam chr2L_dm6 chr2R_dm6 chr3L_dm6 chr3R_dm6 chr4_dm6 chrM_dm6 chrY_dm6 chrX_dm6 > ${i}_dm6_sort.bam

#Extraction of reads count
Cnt_RPE1hap1=$(samtools view -F 0x04 -c ${i}_RPE1v1.1_hap1_sort.bam)
Cnt_RPE1hap2=$(samtools view -F 0x04 -c ${i}_RPE1v1.1_hap2_sort.bam)
Cnt_dm6=$(samtools view -F 0x04 -c ${i}_dm6_sort.bam)

#Creating a file of reads count for each sample
echo -e "Sample\tRPE1hap1\tRPE1hap2\tdm6" > ReadCount/RC_${i}.txt
echo -e "${i}\t${Cnt_RPE1hap1}\t${Cnt_RPE1hap2}\t${Cnt_dm6}" >> ReadCount/RC_${i}.txt

#Indexing
samtools index -@ 8 ${i}_RPE1v1.1_hap1_sort.bam
samtools index -@ 8 ${i}_RPE1v1.1_hap2_sort.bam
samtools index -@ 8 ${i}_dm6_sort.bam

mv ${i}*bam* Bam/.

#Removing variables
unset Cnt_RPE1hap1
unset Cnt_RPE1hap2
unset Cnt_dm6

done