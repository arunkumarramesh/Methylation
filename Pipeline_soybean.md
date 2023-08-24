# Pipeline for getting SMPs and SNPs for Soybean

1. Get data
```
/data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch -X 9999999999999 --option-file  PRJNA432760_soybean_wgbs_SRR_Acc_List.txt  -O data
/data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch -X 9999999999999 --option-file  PRJNA432760_soybean_wgs_SRR_Acc_List.txt  -O data_wgs
for file in *.sra; do /data/proj2/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done

/data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes
grep 'exon' GCF_000004515.6_Glycine_max_v4.0_genomic.gtf | cut -f 1,4,5 >gene_pos.bed
```

2. Trim data
```
for file in *_1.fastq.gz; do java -jar /data/proj2/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 15 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

3. Index genome
```
samtools faidx GCF_000004515.6_Glycine_max_v4.0_genomic.fna
picard CreateSequenceDictionary -R GCF_000004515.6_Glycine_max_v4.0_genomic.fna -O GCF_000004515.6_Glycine_max_v4.0_genomic.dict
```

4. Map genomic data
```
for file in *_1.paired.fq.gz; do bwa mem -t 24 /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna  $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} 2>paired_map_err | samtools view -Sb - > ${file/_1.paired.fq.gz/_mapped.bam} 2>paired_map_err; done
```

5. Add readgroups, filter reads, sort reads, markduplicates and index bams
```
for file in *_mapped.bam ; do picard AddOrReplaceReadGroups -I $file -O ${file/_mapped.bam/_readgroup.bam} -LB species -PL illumina -PU 1 -SM $file; done
for file in *_readgroup.bam; do samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b $file >${file/_readgroup.bam/_mq20.bam}; done
for file in *_mq20.bam; do samtools sort -@ 5 $file >${file/_mq20.bam/_sort.bam} ; done
for file in *_sort.bam ; do picard MarkDuplicates -I $file -O ${file/_sort.bam/_marked.bam} -M ${file/_sort.bam/_metrics.txt} ; done
for file in *_marked.bam ; do picard BuildBamIndex -I $file; done
```

6. Call variants in each sample
```
#!/usr/bin/perl -w


my $dir = "/data/proj2/popgen/a.ramesh/projects/methylomes/soybean/data_wgs";
my $input = "lc_files2_3";
chdir("$dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input") || die "couldn't open the input file $input!";
                    while (my $line = <INPUT_FILE>) {
                        chomp $line;
my @array = split(/\t/,$line);
my $sample = $array[0];
my $file = $array[1];

my $SLURM_header = <<"SLURM";
#\$-cwd
# This option tells gridware to change to the current directory before executing the job
# (default is the home of the user).

#\$-pe serial 2
# Specify this option only for multithreaded jobs that use more than one cpu core.
# IMPORTANT! Don't use more than 4 cores to keep free space for other students!

#\$-l vf=2g
# This option declares the amount of memory the job will use. Specify this value as accurate as possible.
# IMPORTANT! memory request is multiplied by the number of requested CPU cores specified with the -pe.
# Thus, you should divide the overall memory consumption of your job by the number of parallel threads.

#\$-N lc

source ~/.bashrc

cd /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/data_wgs

SLURM

  my $tmp_file = "lc.$sample";
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "gatk HaplotypeCaller -R /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna  -I $file -O  $sample".".g.vcf.gz -ERC GVCF\n";
  close SLURM;
  system("qsub $tmp_file");

}

            close(INPUT_FILE);

```

7. Combine gvcfs
```
gatk CombineGVCFs -R /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna --variant ACC-001.g.vcf.gz --variant ACC-051.g.vcf.gz --variant ACC-076.g.vcf.gz --variant ACC-082.g.vcf.gz --variant ACC-089.g.vcf.gz --variant ACC-100.g.vcf.gz --variant ACC-103.g.vcf.gz --variant ACC-108.g.vcf.gz --variant ACC-1279.g.vcf.gz --variant ACC-128.g.vcf.gz --variant ACC-1297.g.vcf.gz --variant ACC-130.g.vcf.gz --variant ACC-1379.g.vcf.gz --variant ACC-1399.g.vcf.gz --variant ACC-1412.g.vcf.gz --variant ACC-1433.g.vcf.gz --variant ACC-171.g.vcf.gz --variant ACC-172.g.vcf.gz --variant ACC-176.g.vcf.gz --variant ACC-1843.g.vcf.gz --variant ACC-1942.g.vcf.gz --variant ACC-212.g.vcf.gz --variant ACC-215.g.vcf.gz --variant ACC-2210.g.vcf.gz --variant ACC-2218.g.vcf.gz --variant ACC-2219.g.vcf.gz --variant ACC-2225.g.vcf.gz --variant ACC-2226.g.vcf.gz --variant ACC-2227.g.vcf.gz --variant ACC-2228.g.vcf.gz --variant ACC-2229.g.vcf.gz --variant ACC-2230.g.vcf.gz --variant ACC-244.g.vcf.gz --variant ACC-248.g.vcf.gz --variant ACC-250.g.vcf.gz --variant ACC-253.g.vcf.gz --variant ACC-262.g.vcf.gz --variant ACC-281.g.vcf.gz --variant ACC-319.g.vcf.gz --variant ACC-433.g.vcf.gz --variant ACC-438.g.vcf.gz --variant ACC-504.g.vcf.gz --variant ACC-546.g.vcf.gz --variant ACC-616.g.vcf.gz -O soybean.cohort.g.vcf.gz
```

8. genotype samples, select SNP variants and invariant sites
```
gatk --java-options "-Xmx4g" GenotypeGVCFs -R /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna -V soybean.cohort.g.vcf.gz -O soybean.output.vcf.gz

gatk --java-options "-Xmx4g" GenotypeGVCFs -R /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna -all-sites -V soybean.cohort.g.vcf.gz -O soybean.var_invar.vcf.gz
genotypegvcf.sh (END)
```

9. select SNP variants and invariant sites
```
gatk SelectVariants -V soybean.output.vcf.gz -select-type SNP -O soybean.snps.vcf.gz
vcftools --gzvcf soybean.snps.vcf.gz --out soybean_snps_filtered --recode --recode-INFO-all --minDP 20 --minGQ 30  --minQ 30

gatk SelectVariants -V soybean.var_invar.vcf.gz -select-type NO_VARIATION -O soybean.invar.vcf.gz 
vcftools --gzvcf soybean.invar.vcf.gz --out soybean.invar --recode --recode-INFO-all --minDP 20  --minQ 30 --max-missing 0.5 --bed /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/gene_pos.bed
```

10. substitute SNPs into ref genome to create sample specific ref genomes
```
bgzip soybean_snps_filtered.recode.vcf
tabix soybean_snps_filtered.recode.vcf.gz
cat sample_chr | while read -r value1 value2 remainder ; do samtools faidx /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes/GCF_000004515.6_Glycine_max_v4.0_genomic.fna $value2 | bcftools consensus -M N -s $value1 -p ${value1/$/_} -H 1pIu soybean_snps_filtered.recode.vcf.gz >>${value2/$/_soybean.fa}; done
```

11. index each new ref genome and move into seperate folders
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

12. index genome using bismark for methylation mapping
```
/data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark_genome_preparation --hisat2 --verbose --path_to_aligner /data/proj2/popgen/a.ramesh/software/hisat2-2.2.1/ /data/proj2/popgen/a.ramesh/projects/methylomes/soybean/genomes
cat samplenames2 | while read line ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark_genome_preparation --hisat2 --verbose --path_to_aligner /data/proj2/popgen/a.ramesh/software/hisat2-2.2.1/ $line ; done
```

13. May methylomes to sample specific reference genomes
```
sed 's/^/\/data\/proj2\/popgen\/a.ramesh\/projects\/methylomes\/soybean\/pseudogenomes\//' samplenames | paste samplenames - >samplenames4
cat samplenames4 |  while read -r value1 value2 remainder ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /data/proj2/popgen/a.ramesh/software/hisat2-2.2.1/ --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done
```

14. Deduplicate methyalation reads
```
for file in *_bismark_hisat2_pe.bam; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done
```

15. Call methylation variants
```
cat samplenames4 |  while read -r value1 value2 remainder ; do /data/proj2/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam ; done
```

16. Do binomial test
```
setwd("/data/proj2/popgen/a.ramesh/projects/methylomes/soybean/data")
library(dplyr)
library(ggplot2)

covfile <- read.table(file="covfiles")
contextfiles <- read.table(file="contextfiles")

cov <- read.table(file=covfile[1,], header=F)
colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
cov$ID <- paste(cov$chromosome,cov$position,sep = "-")

context <- read.table(file=contextfiles[1,],header=F)
colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
context$ID <- paste(context$chromosome,context$position,sep = "-")

cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
dim(cov_context)
cov_context$chromosome <- gsub("_mapped.bam","", gsub(gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile$V1)[1],"",cov_context$chromosome))
cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]

cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "-")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
cov_context$pval <- 0
noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "NC_007942.1",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "NC_007942.1",]$count.total)

#cov_context <- cov_context[cov_context$count.total > 9,]

b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

cov_context <- cov_context[c(1,2,7,8,9,12)]
cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
cov_context$call <- "U"
cov_context[cov_context$fdr < 0.01,]$call <- "M"
cov_context <- cov_context[-c(6,7)]

samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile[1,])
colnames(cov_context)[ncol(cov_context)] <- samplename
cov_context2 <- cov_context

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles")
  contextfiles <- read.table(file="contextfiles")
  
  cov <- read.table(file=covfile[i,], header=F)
  colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
  cov$ID <- paste(cov$chromosome,cov$position,sep = "-")
  
  context <- read.table(file=contextfiles[i,],header=F)
  colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
  context$ID <- paste(context$chromosome,context$position,sep = "-")
  
  cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
  dim(cov_context)
  cov_context$chromosome <- gsub("_mapped.bam","", gsub(gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile$V1)[i],"",cov_context$chromosome))
  cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]
  
  cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "-")
  cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
  cov_context$pval <- 0
  noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "NC_007942.1",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "NC_007942.1",]$count.total)
  
  #cov_context <- cov_context[cov_context$count.total > 9,]
  
  b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

  cov_context <- cov_context[c(1,2,7,8,9,12)]
  cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
  cov_context$call <- "U"
  cov_context[cov_context$fdr < 0.01,]$call <- "M"
  cov_context <- cov_context[-c(6,7)]
  
  samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile[i,])
  colnames(cov_context)[ncol(cov_context)] <- samplename
  cov_context <- cov_context[c(3,6)]
  cov_context2 <- full_join(cov_context2,cov_context,by="ID")
}

cov_context2$chromosome <- gsub("-.*","",cov_context2$ID)
cov_context2$position   <- as.numeric(gsub(".*-","",cov_context2$ID))
write.table(cov_context2,file="cov_context3.txt",row.names = F)

cov_context3 <- read.table(file="cov_context3.txt",header=T)
na_count <- apply(cov_context3[6:ncol(cov_context3)], 1, function(x) sum(is.na(x)))
na_count <- na_count/ncol(cov_context3[6:ncol(cov_context3)])
cov_context3 <- cov_context3[na_count < 0.5,] # change to appropriate number
cov_context4 <- cov_context3
ploymorphic <- apply(cov_context3[6:ncol(cov_context3)], 1, table)
ploymorphic <- sapply(ploymorphic,length)
cov_context3 <- cov_context3[ploymorphic > 1,]

meta <- cov_context3[1:3]
colnames(meta) <- c("#CHROM","POS","ID")
meta$REF <- "A"
meta$ALT <- "T"
meta$QUAL <- 4000
meta$FILTER <- "PASS"
meta$INFO <- "DP=1000"
meta$FORMAT <- "GT"
cov_context3 <- cov_context3[-c(1:5)]
cov_context3[cov_context3 == "U"] <- "0/0"
cov_context3[cov_context3 == "M"] <- "1/1"
cov_context3[is.na(cov_context3)] <- "./."
write.table(cbind(meta,cov_context3),file="soybean_meth.vcf",quote = F, row.names = F,sep="\t")

meta2 <- cov_context4[1:3]
colnames(meta2) <- c("#CHROM","POS","ID")
meta2$REF <- "A"
meta2$ALT <- "T"
meta2$QUAL <- 4000
meta2$FILTER <- "PASS"
meta2$INFO <- "DP=1000"
meta2$FORMAT <- "GT"
cov_context4 <- cov_context4[-c(1:5)]
cov_context4[cov_context4 == "U"] <- "0/0"
cov_context4[cov_context4 == "M"] <- "1/1"
cov_context4[is.na(cov_context4)] <- "./."
write.table(cbind(meta2,cov_context4),file="soybean_meth_var_invar.vcf",quote = F, row.names = F,sep="\t")
```

17. Some general stats
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
