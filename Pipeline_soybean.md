

```
/data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch -X 9999999999999 --option-file  PRJNA432760_soybean_wgbs_SRR_Acc_List.txt  -O data
/data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch -X 9999999999999 --option-file  PRJNA432760_soybean_wgs_SRR_Acc_List.txt  -O data_wgs
for file in *.sra; do /data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done
```

```
for file in *_1.fastq.gz; do java -jar /data/proj2/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 15 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

```
for file in *_1.paired.fq.gz; do bwa mem -t 24 /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} 2>paired_map_err | samtools view -Sb - > ${file/_1.paired.fq.gz/_mapped.bam} 2>paired_map_err; done
```

```
for file in *_mapped.bam ; do picard AddOrReplaceReadGroups -I $file -O ${file/_mapped.bam/_readgroup.bam} -LB species -PL illumina -PU 1 -SM $file; done
for file in *_readgroup.bam; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b $file >${file/_readgroup.bam/_mq20.bam}; done
for file in *_mq20.bam; do samtools sort -@ 5 $file >${file/_mq20.bam/_sort.bam} ; done
for file in *_sort.bam ; do picard MarkDuplicates -I $file -O ${file/_sort.bam/_marked.bam} -M ${file/_sort.bam/_metrics.txt} ; done
for file in *_marked.bam ; do picard BuildBamIndex -I $file; done
```

```
gatk --java-options "-Xmx4g" GenotypeGVCFs -R /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna -V soybean.cohort.g.vcf.gz -O soybean.output.vcf.gz

gatk --java-options "-Xmx4g" GenotypeGVCFs -R /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna -all-sites -V soybean.cohort.g.vcf.gz -O soybean.var_invar.vcf.gz
genotypegvcf.sh (END)
```

```
gatk SelectVariants -V soybean.output.vcf.gz -select-type SNP -O soybean.snps.vcf.gz
vcftools --gzvcf soybean.snps.vcf.gz --out soybean_snps_filtered --recode --recode-INFO-all --minDP 20 --minGQ 30  --minQ 30

gatk SelectVariants -V soybean.var_invar.vcf.gz -select-type NO_VARIATION -O soybean.invar.vcf.gz 
vcftools --gzvcf soybean.invar.vcf.gz --out soybean.invar --recode --recode-INFO-all --minDP 20  --minQ 30 --max-missing 0.5 --bed /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
```

```
bgzip soybean_snps_filtered.recode.vcf
tabix soybean_snps_filtered.recode.vcf.gz
cat sample_chr | while read -r value1 value2 remainder ; do samtools faidx /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna $value2 | bcftools consensus -M N -s $value1 -p ${value1/$/_} -H 1pIu soybean_snps_filtered.recode.vcf.gz >>${value2/$/_soybean.fa}; done
```

```
cat chrlist | while read line; do /data/proj2/popgen/a.ramesh/software/faSplit byname $line ${line/^/chr_} ; done
cat samplenames | while read line ; do cat $line*.fa >$line.merged.fa  ; done
for file in *.merged.fa ; do samtools faidx $file; done
for file in *.merged.fa ; do picard CreateSequenceDictionary -R $file -O ${file/.fa/.dict}; done
cat samplenames | while read line ; do mkdir $line  ; done
cat samplenames | while read line ; do mv $line.merged.fa $line  ; done
cat samplenames | while read line ; do mv $line.merged.dict $line  ; done
cat samplenames | while read line ; do mv $line.merged.fa.fai $line  ; done
```

```
/data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark_genome_preparation --hisat2 --verbose --path_to_aligner /data/proj2/popgen/a.ramesh/software/hisat2-2.2.1/ /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes
cat samplenames2 | while read line ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark_genome_preparation --hisat2 --verbose --path_to_aligner /data/proj2/popgen/a.ramesh/software/hisat2-2.2.1/ $line ; done
```

```
cat samplenames4 |  while read -r value1 value2 remainder ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /data/proj2/popgen/a.ramesh/software/hisat2-2.2.1/ --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done
```

```
for file in *_bismark_hisat2_pe.bam; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done
```

```
cat samplenames4 |  while read -r value1 value2 remainder ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam ; done
```

```
vcftools --gzvcf soybean_snps_filtered.recode.vcf.gz --out soybean_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --freq --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
vcftools --gzvcf soybean_snps_filtered.recode.vcf.gz --out soybean_snp --min-alleles 2  --max-missing 0.5 --window-pi  10000 --bed /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
vcftools --vcf soybean.invar.recode.vcf --out soybean_invar  --max-missing 0.5 --site-pi --bed /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
vcftools --gzvcf soybean_snps_filtered.recode.vcf.gz --out soybean_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --het --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
vcftools --gzvcf soybean_snps_filtered.recode.vcf.gz --out soybean_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --geno-r2 --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed --ld-window-bp 100

vcftools --gzvcf soybean_snps_filtered.recode.vcf.gz --out soybean_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --recode --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed

cd /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/data/
vcftools --vcf soybean_meth_all.vcf  --out soybean_meth --min-alleles 2 --max-alleles 2 --max-missing 0.5 --freq --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
vcftools --vcf soybean_meth_all.vcf  --out soybean_meth --min-alleles 2 --max-alleles 2 --max-missing 0.5 --site-pi --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
vcftools --vcf soybean_meth_all.vcf  --out soybean_meth --min-alleles 2 --max-alleles 2 --max-missing 0.5 --geno-r2 --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed --ld-window-bp 100
```
