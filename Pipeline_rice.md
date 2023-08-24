# Pipeline for getting SMPs and SNPs for rice

1. Download data from NCBI
```
## Genomes usually donwloaded from Ensembl, need to ensure it includes the choloplast contig. Index with samtools faidx.

cd /proj/popgen/a.ramesh/projects/methylomes/rice
/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file PRJNA597475_rice_Acc_List.txt  -O data
/proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file SRR_Acc_List_rice_rnaseq.txt  -O data_rna

cd /proj/popgen/a.ramesh/projects/methylomes/rice/data
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done

cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna
for file in *.sra; do /proj/popgen/a.ramesh/software/sratoolkit.3.0.0-centos_linux64/bin/fastq-dump --split-3 --gzip  $file; done

```

2. Trim data
```
cd /proj/popgen/a.ramesh/projects/methylomes/lyrata/data
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done

cd /proj/popgen/a.ramesh/projects/methylomes/lyrata/data_rna
for file in *_1.fastq.gz; do java -jar /proj/popgen/a.ramesh/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 $file ${file/_1.fastq.gz/_2.fastq.gz} ${file/_1.fastq.gz/_1.paired.fq.gz} ${file/_1.fastq.gz/_1.unpaired.fq.gz} ${file/_1.fastq.gz/_2.paired.fq.gz} ${file/_1.fastq.gz/_2.unpaired.fq.gz} ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36; done
```

3. Index reference genome and map RNA-seq data with with STAR 2-pass mode

```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/genomes
/proj/popgen/a.ramesh/software/STAR-2.7.10b/source/STAR --runMode genomeGenerate --genomeDir STARindex/ --genomeFastaFiles /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --sjdbGTFfile Oryza_sativa.IRGSP-1.0.55.gtf --sjdbOverhang 149 --runThreadN 20

cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna/
ulimit -n 10000
for file in *_1.paired.fq.gz; do /proj/popgen/a.ramesh/software/STAR-2.7.10b/source/STAR --runMode alignReads --genomeDir /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/STARindex/ --readFilesIn $file ${file/_1.paired.fq.gz/_2.paired.fq.gz} --readFilesCommand zcat --outFileNamePrefix ${file/_1.paired.fq.gz//}  --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --twopassMode Basic --runThreadN 20  ; done

find . -type f -name "Aligned.sortedByCoord.out.bam" -printf "/%P\n" | while read FILE ; do DIR=$(dirname "$FILE" ); mv ."$FILE" ."$DIR""$DIR".sort.bam;done
find . -name '*.sort.bam' -exec mv {} . \;

```

4. Add read groupsls
```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna/

for file in SRR*.sort.bam ; do java -jar /proj/popgen/a.ramesh/software/picard.jar AddOrReplaceReadGroups -I $file -O ${file/.sort.bam/_readgroup.bam} -LB species -PL illumina -PU 1 -SM $file; done

```

5. Mark duplicate reads and index
```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna/

for file in SRR*.sort.bam ; do java -jar /proj/popgen/a.ramesh/software/picard.jar AddOrReplaceReadGroups -I $file -O ${file/.sort.bam/_readgroup.bam} -LB species -PL illumina -PU 1 -SM $file; done

java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751892_readgroup.bam -I SRR10751893_readgroup.bam -O Minghui63_marked.bam -M Minghui63_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751900_readgroup.bam -I SRR10751901_readgroup.bam -O ZhenShan97_marked.bam -M ZhenShan97_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751908_readgroup.bam -I SRR10751909_readgroup.bam -O Nipponare_marked.bam -M Nipponare_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751916_readgroup.bam -I SRR10751917_readgroup.bam -O HKG98_marked.bam -M HKG98_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751918_readgroup.bam -I SRR10751919_readgroup.bam -O Garia_marked.bam -M Garia_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751920_readgroup.bam -I SRR10751921_readgroup.bam -O Aijiaonante_marked.bam -M Aijiaonante_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751922_readgroup.bam -I SRR10751923_readgroup.bam -O Xibaizhan_marked.bam -M Xibaizhan_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751924_readgroup.bam -I SRR10751925_readgroup.bam -O WH139_marked.bam -M WH139_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751926_readgroup.bam -I SRR10751927_readgroup.bam -O Nanjing11_marked.bam -M Nanjing11_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751928_readgroup.bam -I SRR10751929_readgroup.bam -O 9311_marked.bam -M 9311_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751930_readgroup.bam -I SRR10751931_readgroup.bam -O Y134_marked.bam -M Y134_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751932_readgroup.bam -I SRR10751933_readgroup.bam -O IR72_marked.bam -M IR72_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751934_readgroup.bam -I SRR10751935_readgroup.bam -O PeiC122_marked.bam -M PeiC122_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751936_readgroup.bam -I SRR10751937_readgroup.bam -O YongChalByo_marked.bam -M YongChalByo_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751938_readgroup.bam -I SRR10751939_readgroup.bam -O TAINO38_marked.bam -M TAINO38_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751940_readgroup.bam -I SRR10751941_readgroup.bam -O Ginga_marked.bam -M Ginga_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751942_readgroup.bam -I SRR10751943_readgroup.bam -O LABELLE_marked.bam -M LABELLE_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751944_readgroup.bam -I SRR10751945_readgroup.bam -O Latisai1_marked.bam -M Latisai1_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751946_readgroup.bam -I SRR10751947_readgroup.bam -O GasymHany_marked.bam -M GasymHany_metrics.txt
java -jar /proj/popgen/a.ramesh/software/picard.jar MarkDuplicates -I SRR10751948_readgroup.bam -I SRR10751949_readgroup.bam -O BASMATI385_marked.bam -M BASMATI385_metrics.txt

for file in *_marked.bam ; do java -jar /proj/popgen/a.ramesh/software/picard.jar BuildBamIndex -I $file ; done

```

6. SplitNCigarReads
```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna/
for file in *_marked.bam ; do /proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk SplitNCigarReads -R /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa  -I $file -O  ${file/_marked.bam/_split.bam}  ; done
for file in *_split.bam ; do picard BuildBamIndex -I $file; done
ls *_split.bam > lc_files
ls *_split.bam | sed 's/*_split.bam//' | paste - lc_files >lc_files2 ## for haplotyper caller script (file also called hc_files2)
```

7. Run haplotype caller using perl script to parallelize
```
#!/usr/bin/perl -w


my $dir = "/proj/popgen/a.ramesh/projects/methylomes/rice/data_rna";
my $input = "lc_files2";
chdir("$dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input") || die "couldn't open the input file $input!";
                    while (my $line = <INPUT_FILE>) {
                        chomp $line;
my @array = split(/\t/,$line);
my $sample = $array[0];
my $file = $array[1];

my $SLURM_header = <<"SLURM";
#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=100000
#SBATCH -J lc
#SBATCH -o lc.out
#SBATCH -e lc.err

cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna

SLURM

  my $tmp_file = "hc.$sample";
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk HaplotypeCaller -R /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa  -I $file -O  $sample".".g.vcf.gz  -ERC GVCF\n";
  close SLURM;
  system("sbatch $tmp_file");

}

            close(INPUT_FILE);

```


8. Combine GVCFs
```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna/
/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk CombineGVCFs -R /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --variant C019.g.vcf.gz --variant C051.g.vcf.gz --variant C135.g.vcf.gz --variant C139.g.vcf.gz --variant C148.g.vcf.gz --variant C151.g.vcf.gz --variant MH63.g.vcf.gz --variant NIP.g.vcf.gz --variant W081.g.vcf.gz --variant W105.g.vcf.gz --variant W125.g.vcf.gz --variant W128.g.vcf.gz --variant W161.g.vcf.gz --variant W169.g.vcf.gz --variant W257.g.vcf.gz --variant W261.g.vcf.gz --variant W286.g.vcf.gz --variant W294.g.vcf.gz --variant W306.g.vcf.gz --variant ZS97.g.vcf.gz -O rice.cohort.g.vcf.gz
```

9. Genotype combined GVCFs

```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna/
/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk --java-options "-Xmx4g" GenotypeGVCFs -R /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa  -V rice.cohort.g.vcf.gz -O rice.output.vcf.gz

/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk --java-options "-Xmx4g" GenotypeGVCFs -R /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -all-sites -V rice.cohort.g.vcf.gz -O rice.var_invar.vcf.gz
```

10. Filter vcf, only keep high quality SNPs
```
#/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk SelectVariants -V rice.output.vcf.gz -select-type SNP -O rice.snps.vcf.gz
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf rice.snps.vcf.gz --out rice_snps_filtered --recode --recode-INFO-all --minDP 20 --minGQ 30 --minQ 30
#sed -i -e 's/9311_split.bam/C148/' -e 's/Aijiaonante_split.bam/C019/' -e 's/BASMATI385_split.bam/W306/' -e 's/Garia_split.bam/W286/' -e 's/GasymHany_split.bam/W081/' -e 's/Ginga_split.bam/W261/' -e 's/HKG98_split.bam/W105/' -e 's/IR72_split.bam/W169/' -e 's/LABELLE_split.bam/W294/' -e 's/Latisai1_split.bam/W257/' -e 's/Minghui63_split.bam/MH63/' -e 's/Nanjing11_split.bam/C151/' -e 's/Nipponare_split.bam/NIP/' -e 's/PeiC122_split.bam/C051/' -e 's/TAINO38_split.bam/W128/' -e 's/WH139_split.bam/C139/' -e 's/Xibaizhan_split.bam/C135/' -e 's/Y134_split.bam/W161/' -e 's/YongChalByo_split.bam/W125/' -e 's/ZhenShan97_split.bam/ZS97/' rice_snps_filtered.recode.vcf


#/proj/popgen/a.ramesh/software/gatk-4.3.0.0/gatk SelectVariants -V rice.var_invar.vcf.gz -select-type NO_VARIATION -O rice.invar.vcf.gz
#/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf rice.invar.vcf.gz --out rice.invar --recode --recode-INFO-all --minDP 20 --minGQ 30 --minQ 30 --max-missing 0.5 --bed /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/gene_pos.bed
#sed -i -e 's/9311_split.bam/C148/' -e 's/Aijiaonante_split.bam/C019/' -e 's/BASMATI385_split.bam/W306/' -e 's/Garia_split.bam/W286/' -e 's/GasymHany_split.bam/W081/' -e 's/Ginga_split.bam/W261/' -e 's/HKG98_split.bam/W105/' -e 's/IR72_split.bam/W169/' -e 's/LABELLE_split.bam/W294/' -e 's/Latisai1_split.bam/W257/' -e 's/Minghui63_split.bam/MH63/' -e 's/Nanjing11_split.bam/C151/' -e 's/Nipponare_split.bam/NIP/' -e 's/PeiC122_split.bam/C051/' -e 's/TAINO38_split.bam/W128/' -e 's/WH139_split.bam/C139/' -e 's/Xibaizhan_split.bam/C135/' -e 's/Y134_split.bam/W161/' -e 's/YongChalByo_split.bam/W125/' -e 's/ZhenShan97_split.bam/ZS97/' rice.invar.recode.vcf

/proj/popgen/a.ramesh/software/vcftools-vcftools-581c231/bin/vcftools --gzvcf rice.var_invar.vcf.gz --out rice.allsites --recode --recode-INFO-all --minDP 20 --minQ 30 --max-missing 1 --bed /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/gene_pos.bed
sed -i -e 's/9311_split.bam/C148/' -e 's/Aijiaonante_split.bam/C019/' -e 's/BASMATI385_split.bam/W306/' -e 's/Garia_split.bam/W286/' -e 's/GasymHany_split.bam/W081/' -e 's/Ginga_split.bam/W261/' -e 's/HKG98_split.bam/W105/' -e 's/IR72_split.bam/W169/' -e 's/LABELLE_split.bam/W294/' -e 's/Latisai1_split.bam/W257/' -e 's/Minghui63_split.bam/MH63/' -e 's/Nanjing11_split.bam/C151/' -e 's/Nipponare_split.bam/NIP/' -e 's/PeiC122_split.bam/C051/' -e 's/TAINO38_split.bam/W128/' -e 's/WH139_split.bam/C139/' -e 's/Xibaizhan_split.bam/C135/' -e 's/Y134_split.bam/W161/' -e 's/YongChalByo_split.bam/W125/' -e 's/ZhenShan97_split.bam/ZS97/' rice.allsites.recode.vcf
```

11. Get chromosome and sample names, can be done many ways
```
cut -f 1 ../genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.fai >chrlist
ls *_split.bam | sed 's/_split.bam//' >samplenames
```

12. Create combined list containing every pair of chromosome and sample names. In R. 
```
file1 <- read.table(file="samplenames")
file2 <- read.table(file="chrlist")

fileall <- expand.grid(file1$V1, file2$V1)  
write.table(fileall, file = "sample_chr", quote = F, sep = "\t", row.names = F, col.names = F)
```

13. Create pseudogenomes by substituting SNPs for each sample from VCF into the reference genome
```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data_rna/
/proj/popgen/a.ramesh/software/htslib-1.16/bgzip rice_snps_filtered.recode.vcf
/proj/popgen/a.ramesh/software/htslib-1.16/tabix  rice_snps_filtered.recode.vcf.gz
cat sample_chr | while read -r value1 value2 remainder ; do samtools faidx /proj/popgen/a.ramesh/projects/methylomes/rice/genomes/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa $value2 | /proj/popgen/a.ramesh/software/bcftools-1.16/bcftools consensus -M N -s $value1 -p ${value1/$/_} -H 1pIu rice_snps_filtered.recode.vcf.gz >>${value2/$/_rice.fa}; done
cat chrlist | while read line; do /data/proj2/popgen/a.ramesh/software/faSplit byname $line ${line/^/chr_} ; done
mkdir /proj/popgen/a.ramesh/projects/methylomes/rice/pseudogenomes
mv *fa /proj/popgen/a.ramesh/projects/methylomes/rice/pseudogenomes
cd /proj/popgen/a.ramesh/projects/methylomes/rice/pseudogenomes
cat samplenames | while read line ; do cat $line*.fa >$line.merged.fa  ; done
for file in *.merged.fa ; do samtools faidx $file; done
for file in *.merged.fa ; do java -jar /proj/popgen/a.ramesh/software/picard.jar CreateSequenceDictionary -R $file -O ${file/.fa/.dict}; done
cat samplenames | while read line ; do mkdir $line  ; done
cat samplenames | while read line ; do mv $line.merged.fa $line  ; done
cat samplenames | while read line ; do mv $line.merged.dict $line  ; done
cat samplenames | while read line ; do mv $line.merged.fa.fai $line  ; done
cat samplenames2 | while read line ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_genome_preparation --hisat2 --verbose --path_to_aligner /proj/popgen/a.ramesh/software/hisat2-2.2.1/ $line ; done

```

14. Now map methylation reads with Bismark
```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data

sed 's/^/\/proj\/popgen\/a.ramesh\/projects\/methylomes\/rice\/pseudogenomes\//' samplenames | paste samplenames >samplenames2

cat samplenames2 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark --multicore 4 --hisat2 --path_to_hisat2 /proj/popgen/a.ramesh/software/hisat2-2.2.1/ --genome_folder $value2 -1 $value1.1.paired.fq.gz -2 $value1.2.paired.fq.gz  ; done

```

15. Deduplicate methylation reads with Bismark
```
cd /proj/popgen/a.ramesh/projects/methylomes/rice/data
for file in *.1.paired_bismark_hisat2_pe.bam ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/deduplicate_bismark --bam $file ; done
```

16. Call methylation variants
```
d /proj/popgen/a.ramesh/projects/methylomes/rice/data
cat samplenames2 |  while read -r value1 value2 remainder ; do /proj/popgen/a.ramesh/software/Bismark-0.24.0/bismark_methylation_extractor --multicore 4 --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder $value2 $value1.1.paired_bismark_hisat2_pe.deduplicated.bam ; done
```

17. This script is to create methylation vcfs. I need add this after annotating it properly. Lots done here. In R.
```
setwd("/data/proj2/popgen/a.ramesh/projects/methylomes/rice/data")
library(dplyr)
library(ggplot2)

covfile <- read.table(file="covfiles")
contextfiles <- read.table(file="contextfiles")

cov <- read.table(file=covfile[1,], header=F)
colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
cov$ID <- paste(cov$chromosome,cov$position,sep = "_")

context <- read.table(file=contextfiles[1,],header=F)
colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
context$ID <- paste(context$chromosome,context$position,sep = "_")

cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
dim(cov_context)
cov_context$chromosome <- gsub(gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile$V1)[1],"",cov_context$chromosome)
cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]

cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "_")
cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
cov_context$pval <- 0
noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "Pt",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "Pt",]$count.total)
if (is.na(noncoversionrate)) {
  noncoversionrate <- 0
}
samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile[1,])

b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

cov_context <- cov_context[c(1,2,7,8,9,12)]
cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
cov_context$call <- "U"
cov_context[cov_context$fdr < 0.01,]$call <- "M"
cov_context <- cov_context[-c(6,7)]

colnames(cov_context)[ncol(cov_context)] <- samplename
cov_context2 <- cov_context

for (i in 2:nrow(covfile)){
  covfile <- read.table(file="covfiles")
  contextfiles <- read.table(file="contextfiles")
  
  print(i)
  
  cov <- read.table(file=covfile[i,], header=F)
  colnames(cov) <- c("chromosome", "position", "end.position", "methylation.percentage", "count.methylated", "count.unmethylated" )
  cov$ID <- paste(cov$chromosome,cov$position,sep = "_")
  
  context <- read.table(file=contextfiles[i,],header=F)
  colnames(context) <- c("chromosome", "position", "strand", "count.methylated", "count.unmethylated", "C-context", "trinucleotide context")
  context$ID <- paste(context$chromosome,context$position,sep = "_")
  
  cov_context <- inner_join(cov,context[c(3,6:8)],by="ID")
  dim(cov_context)
  cov_context$chromosome <- gsub(gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile$V1)[i],"",cov_context$chromosome)
  cov_context <- cov_context[cov_context$count.methylated + cov_context$count.unmethylated > 4, ]
  
  cov_context$ID <- paste(cov_context$chromosome,cov_context$position,sep = "_")
  cov_context$count.total <- cov_context$count.methylated + cov_context$count.unmethylated
  cov_context$pval <- 0
  noncoversionrate <- sum(cov_context[cov_context$chromosome %in% "Pt",]$count.methylated)/sum(cov_context[cov_context$chromosome %in% "Pt",]$count.total)
  if (is.na(noncoversionrate)) {
    noncoversionrate <- 0
  }
  samplename <- gsub(".1.paired_bismark_hisat2_pe.deduplicated.bismark.cov","",covfile[i,])
    
  b <- apply(cov_context[c(5,11)],1,binom.test,p = noncoversionrate, alternative = c("greater"))
  cov_context$pval <- do.call(rbind,lapply(b,function(v){v$p.value}))

  cov_context <- cov_context[c(1,2,7,8,9,12)]
  cov_context$fdr <- p.adjust(cov_context$pval,method = "fdr")
  cov_context$call <- "U"
  cov_context[cov_context$fdr < 0.01,]$call <- "M"
  cov_context <- cov_context[-c(6,7)]
  

  colnames(cov_context)[ncol(cov_context)] <- samplename
  cov_context <- cov_context[c(3,6)]
  cov_context2 <- full_join(cov_context2,cov_context,by="ID")
}
cov_context2$chromosome <- gsub("_.*","",cov_context2$ID)
cov_context2$position   <- as.numeric(gsub(".*_","",cov_context2$ID))

write.table(cov_context2,file="cov_context3.txt",row.names = F)

etwd("/data/proj2/popgen/a.ramesh/projects/methylomes/rice/data")

#cov_context3 <- read.table(file="cov_context3.txt",header=T)
#na_count <- apply(cov_context3[6:ncol(cov_context3)], 1, function(x) sum(is.na(x)))
#na_count <- na_count/ncol(cov_context3[6:ncol(cov_context3)])
#cov_context3 <- cov_context3[na_count < 0.5,] # change to appropriate number
#cov_context4 <- cov_context3
#ploymorphic <- apply(cov_context3[6:ncol(cov_context3)], 1, table)
#ploymorphic <- sapply(ploymorphic,length)
#cov_context3 <- cov_context3[ploymorphic > 1,]

#meta <- cov_context3[1:3]
#colnames(meta) <- c("#CHROM","POS","ID")
#meta$REF <- "A"
#meta$ALT <- "T"
#meta$QUAL <- 4000
#meta$FILTER <- "PASS"
#meta$INFO <- "DP=1000"
#meta$FORMAT <- "GT"

#cov_context3 <- cov_context3[-c(1:5)]
#cov_context3[cov_context3 == "U"] <- "0/0"
#cov_context3[cov_context3 == "M"] <- "1/1"
#cov_context3[is.na(cov_context3)] <- "./."

#meta <- cbind(meta,cov_context3)
#meta <- meta[meta$`#CHROM` %in% c("1","2","3","4","5","6","7","8","9","10","11","12"),]
#meta$`#CHROM` <- as.numeric(meta$`#CHROM`)
#meta$POS <- as.numeric(meta$POS)
#meta <- meta[order(meta$`#CHROM`,meta$POS),]

#write.table(meta,file="rice_meth.vcf",quote = F, row.names = F,sep="\t")

#meta2 <- cov_context4[1:3]
#colnames(meta2) <- c("#CHROM","POS","ID")
#meta2$REF <- "A"
#meta2$ALT <- "T"
#meta2$QUAL <- 4000
#meta2$FILTER <- "PASS"
#meta2$INFO <- "DP=1000"
#meta2$FORMAT <- "GT"
#cov_context4 <- cov_context4[-c(1:5)]
#cov_context4[cov_context4 == "U"] <- "0/0"
#cov_context4[cov_context4 == "M"] <- "1/1"
#cov_context4[is.na(cov_context4)] <- "./."
#meta2 <- cbind(meta2,cov_context4)
#meta2 <- meta2[meta2$`#CHROM` %in% c("1","2","3","4","5","6","7","8","9","10","11","12"),]
#meta2$`#CHROM` <- as.numeric(meta2$`#CHROM`)
#meta2$POS <- as.numeric(meta2$POS)
#meta2 <- meta2[order(meta2$`#CHROM`,meta2$POS),]
#write.table(meta2,file="rice_meth_var_invar.vcf",quote = F, row.names = F,sep="\t")

#first run 'vcftools --vcf rice_meth_var_invar_all.vcf  --out rice_meth_var_invar --max-missing 0.5  --recode --bed  gene_pos.bed '
vcfnames <- as.character(t(read.table(file="vcfnames")))

cov_context3 <- read.table(file="rice_meth_var_invar.recode.vcf")
colnames(cov_context3) <- vcfnames

meta <- cov_context3[c(1,2,6)]
colnames(meta) <- c("CHROM","POS","VAULE")
cov_context3 <- cov_context3[-c(1:9)]

cov_context3 <- cov_context3[grep("^C",colnames(cov_context3))]

cov_context3[cov_context3 == "./."] <- "CC"
cov_context3[cov_context3 == "0/0"] <- "DD"
cov_context3[cov_context3 == "1/1"] <- "MM"
library(tidyr)
cov_context3 <- unite(cov_context3, col='combined', colnames(cov_context3), sep='')
meta <- cbind(meta,cov_context3)

meta_1 <- meta[meta$CHROM %in% 1,]
meta_2 <- meta[meta$CHROM %in% 2,]
meta_3 <- meta[meta$CHROM %in% 3,]
meta_4 <- meta[meta$CHROM %in% 4,]
meta_5 <- meta[meta$CHROM %in% 5,]
meta_6 <- meta[meta$CHROM %in% 6,]
meta_7 <- meta[meta$CHROM %in% 7,]
meta_8 <- meta[meta$CHROM %in% 8,]
meta_9 <- meta[meta$CHROM %in% 9,]
meta_10 <- meta[meta$CHROM %in% 10,]
meta_11 <- meta[meta$CHROM %in% 11,]
meta_12 <- meta[meta$CHROM %in% 12,]

write.table(meta_1,file="rice_multihetsep_meth_1",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_2,file="rice_multihetsep_meth_2",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_3,file="rice_multihetsep_meth_3",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_4,file="rice_multihetsep_meth_4",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_5,file="rice_multihetsep_meth_5",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_6,file="rice_multihetsep_meth_6",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_7,file="rice_multihetsep_meth_7",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_8,file="rice_multihetsep_meth_8",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_9,file="rice_multihetsep_meth_9",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_10,file="rice_multihetsep_meth_10",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_11,file="rice_multihetsep_meth_11",quote = F, col.names = F,row.names = F,sep="\t")
write.table(meta_12,file="rice_multihetsep_meth_12",quote = F, col.names = F,row.names = F,sep="\t")


```

18. Get vcf header
```
zcat lyrata_snps_filtered.recode.vcf.gz | grep '##' > vcfheader
cat vcfheader lyrata_meth.vcf >lyrata_meth_all.vcf

```

19. Get summary statistics for methylation and genomic variants. Done on biallelic SNPs. Further NA filtering in R. Only keep gene body variants.
```
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/
awk '$3 ~ /^gene$/' Arabidopsis_lyrata.v.1.0.55.chr.gff3 | cut -f 1,4,5 >gene_pos.bed
 
cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data_rna/
vcftools --gzvcf lyrata.snps_filtered.lc.recode.vcf.gz --out lyrata_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --freq --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed
vcftools --gzvcf lyrata.snps_filtered.lc.recode.vcf.gz --out lyrata_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --site-pi --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed
vcftools --gzvcf lyrata.snps_filtered.lc.recode.vcf.gz --out lyrata_snp --min-alleles 2 --max-alleles 2 --max-missing 0.5 --het --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed

cd /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/data/
vcftools --vcf lyrata_meth_all.vcf  --out lyrata_meth --min-alleles 2 --max-alleles 2 --max-missing 0.5 --freq --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed
vcftools --vcf lyrata_meth_all.vcf  --out lyrata_meth --min-alleles 2 --max-alleles 2 --max-missing 0.5 --site-pi --bed  /data/proj2/popgen/a.ramesh/projects/methylomes/lyrata/genomes/gene_pos.bed

```
