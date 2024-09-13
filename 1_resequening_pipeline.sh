#!/bin/bash

ref=/database/Reference/Osmia_excavata/Osmia_excavata.fa
data=/project/home/upload/Public/bi_populaton
/software/gatk-4.1.4.0/gatk CreateSequenceDictionary -R ${ref} -O ${ref}.dict
mkdir bam
#/software/gatk-4.1.4.0/gatk CreateSequenceDictionary -R ${ref} -O ${ref}.dict
#/software/bwa/current/bwa index ${ref}
#samtools faidx ${ref}
for p in *_R1.fq.gz
do
n=`echo $p | sed 's/_R1.fq.gz//g'`

/software/bwa/current/bwa mem -t 100 -R '@RG\tID:$n\tPL:illumina\tSM:$n' ${ref} ${data}/${n}_R1.fq.gz ${data}/${n}_R2.fq.gz >${data}/bam/${n}.sam && echo "${n} mapping is done" #change fq 

samtools view -Sb ${data}/bam/${n}.sam > ${data}/bam/${n}.raw.bam && echo "${n} samt2bam is done"

rm -rf ${data}/bam/${n}.sam && echo "${n} rm_sam is done"

samtools sort -@ 4 -m 4G -O bam -T sorted -o ${data}/bam/${n}.sort.bam ${data}/bam/${n}.raw.bam && echo "/${n}.sorted.bam is done"

/software/gatk-4.1.4.0/gatk MarkDuplicates -I ${data}/bam/${n}.sort.bam -O ${data}/bam/${n}.sorted.markdup.bam -M ${data}/bam/${n}.sorted.markdup_metrics.txt

samtools index ${data}/bam/${n}.sorted.markdup.bam

/software/gatk-4.1.4.0/gatk HaplotypeCaller -R ${ref} -I ${data}/bam/${n}.sorted.markdup.bam  --emit-ref-confidence GVCF -O ${data}/bam/${n}.g.raw.vcf

/software/gatk-4.1.4.0/gatk GenotypeGVCFs -R ${ref} -V ${data}/bam/${n}.g.raw.vcf -O ${data}/bam/${n}.raw.vcf

/software/gatk-4.1.4.0/gatk SelectVariants -V ${data}/bam/${n}.raw.vcf -select-type SNP -O ${data}/bam/${n}.select.snp.vcf

/software/gatk-4.1.4.0/gatk SelectVariants -V ${data}/bam/${n}.raw.vcf -select-type INDEL -O ${data}/bam/${n}.select.indel.vcf

/gatk-4.1.4.0/gatk VariantFiltration -V ${data}/bam/${n}.select.snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "PASS" -O ${data}/bam/${n}.filter.select.snp.vcf

done
