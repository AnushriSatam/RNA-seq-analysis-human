#!/bin/bash

set -euo pipefail

WORKDIR="anushri/rna_seq_analysis"

STAR="anushri/softwares/STAR-2.7.10b/bin/Linux_x86_64/STAR"
STAR_INDEX="${WORKDIR}/results/star_index"
SORTMERNA="softwares/sortmerna_4.3.7/bin/sortmerna"
SORTMERNA_DB_DIR="rna-seq/software/silva"
GTF="anushri/human/gencode.v49.basic.annotation.gtf"
GENOME="anushri/human/GRCh38.primary_assembly.genome.fa"

mkdir -p ${WORKDIR}/raw_fastqc
mkdir -p ${WORKDIR}/results/sortmerna
mkdir -p ${WORKDIR}/results/star
mkdir -p ${WORKDIR}/results/star_index
mkdir -p ${WORKDIR}/results/trimmed_fastqc


# FastQC (Quality check)
for f in ${WORKDIR}/*.fastq.gz
do
    fastqc "$f" -o ${WORKDIR}/raw_fastqc -t 32
done

# fastp (trimming)
mkdir -p ${WORKDIR}/results/fastp

for r1 in ${WORKDIR}/*_R1_001.fastq.gz
do
    BASENAME=$(basename $r1 _R1_001.fastq.gz)

    r2=${WORKDIR}/${BASENAME}_R2_001.fastq.gz

    fastp \
        -i $r1 \
        -I $r2 \
        -o ${WORKDIR}/results/fastp/${BASENAME}_R1_trimmed.fastq.gz \
        -O ${WORKDIR}/results/fastp/${BASENAME}_R2_trimmed.fastq.gz \
        -h ${WORKDIR}/results/fastp/${BASENAME}_fastp.html \
        -j ${WORKDIR}/results/fastp/${BASENAME}_fastp.json \
        -q 30 \
        -l 50 \
        --detect_adapter_for_pe \
        -w 16
done

# Quality check after trimming

for f in ${WORKDIR}/results/fastp/*.fastq.gz
do
    fastqc "$f" -o ${WORKDIR}/results/trimmed_fastqc -t 32
done

# kraken-contamination screening

kraken_db="path_to/Kraken_PlusPFP_16_2025"
kraken_outdir="${WORKDIR}/results/kraken"
mkdir -p ${kraken_outdir}

for r1 in ${WORKDIR}/results/fastp/*_R1_trimmed.fastq.gz
do
    BASENAME=$(basename $r1 _R1_trimmed.fastq.gz)
    r2=${WORKDIR}/results/fastp/${BASENAME}_R2_trimmed.fastq.gz

    kraken2 \
        --db $kraken_db \
        --paired \
        --gzip-compressed \
        --threads 32 \
        --output ${kraken_outdir}/${BASENAME}_kraken2_output.txt \
        --report ${kraken_outdir}/${BASENAME}_kraken2_report.txt \
        --report-minimizer-data \
        $r1 $r2

    echo "Kraken2 done: ${BASENAME}"
done

# SortMeRNA — rRNA removal

for r1 in ${WORKDIR}/results/fastp/*_R1_trimmed.fastq.gz
do
    BASENAME=$(basename $r1 _R1_trimmed.fastq.gz)
    r2=${WORKDIR}/results/fastp/${BASENAME}_R2_trimmed.fastq.gz

    rm -rf /path_to/sortmerna_index/idx/kvdb/* 2>/dev/null || true
    rm -rf /path_to/sortmerna_index/idx/readb/* 2>/dev/null || true

    $SORTMERNA \
        --ref ${SORTMERNA_DB_DIR}/rfam-5s-database-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/rfam-5.8s-database-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-arc-16s-id95.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-arc-23s-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-bac-16s-id90.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-bac-23s-id98.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-euk-18s-id95.fasta \
        --ref ${SORTMERNA_DB_DIR}/silva-euk-28s-id98.fasta \
        --reads $r1 \
        --reads $r2 \
        --aligned ${WORKDIR}/results/sortmerna/${BASENAME}_rRNA \
        --other ${WORKDIR}/results/sortmerna/${BASENAME}_mRNA \
        --fastx \
        --paired-in \
        --out2 \
        --threads 32 \
        --index 0 \
        --workdir /path_to/sortmerna_index/idx

    mv ${WORKDIR}/results/sortmerna/${BASENAME}_mRNA_fwd.fq.gz \
       ${WORKDIR}/results/sortmerna/${BASENAME}_R1_mRNA.fastq.gz

    mv ${WORKDIR}/results/sortmerna/${BASENAME}_mRNA_rev.fq.gz \
       ${WORKDIR}/results/sortmerna/${BASENAME}_R2_mRNA.fastq.gz

    echo "SortMeRNA done: ${BASENAME}"
done

#STAR indexing
#was done only once
${STAR} \
    --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir ${STAR_INDEX} \
    --genomeFastaFiles ${GENOME} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang 150 \
    --genomeSAindexNbases 14

#echo "STAR index built."

#STAR alignment

for r1 in ${WORKDIR}/results/sortmerna/*_R1_mRNA.fastq.gz
do
    BASENAME=$(basename $r1 _R1_mRNA.fastq.gz)
    r2=${WORKDIR}/results/sortmerna/${BASENAME}_R2_mRNA.fastq.gz

    ${STAR} \
        --runThreadN 32 \
        --genomeDir ${STAR_INDEX} \
        --readFilesIn $r1 $r2 \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI NM MD \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 2 \
        --alignIntronMin 20 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        --outFileNamePrefix ${WORKDIR}/results/star/${BASENAME}_ \
        --outBAMsortingThreadN 16

samtools index ${WORKDIR}/results/star/${BASENAME}_Aligned.sortedByCoord.out.bam

    #echo "STAR done: ${BASENAME}"
done


#SALMON transcript counts

SALMON="path_to/salmon/bin/salmon"
SALMON_INDEX="${WORKDIR}/results/salmon_index"
SALMON_OUTDIR="${WORKDIR}/results/salmon"

TRANSCRIPTS_FA="anushri/human/gencode.v49.transcripts.fa.gz"

mkdir -p ${SALMON_INDEX}
mkdir -p ${SALMON_OUTDIR}

#indexing-must be done only once

${SALMON} index \
    -t ${TRANSCRIPTS_FA} \
    -i ${SALMON_INDEX} \
    -k 31
echo "Salmon index built."

#quantification

for r1 in ${WORKDIR}/results/sortmerna/*_R1_mRNA.fastq.gz
do
    BASENAME=$(basename $r1 _R1_mRNA.fastq.gz)
    r2=${WORKDIR}/results/sortmerna/${BASENAME}_R2_mRNA.fastq.gz

    ${SALMON} quant \
        -i ${SALMON_INDEX} \
        -l A \
        -1 $r1 \
        -2 $r2 \
        -p 16 \
        --validateMappings \
        -o ${SALMON_OUTDIR}/${BASENAME}_quant

    echo "Salmon done: ${BASENAME}"
done
