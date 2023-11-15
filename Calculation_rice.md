# 1. extract Indica and filter
## 1.1 SNP
####  extract indica
	cd new_dataset/pops
	vcftools --vcf /data/home/students/b.huo/new_dataset/filtering_steps/snps_chr.vcf --keep Indica1.txt --max-missing 0.8 --maf 0.2 --recode --out indica1
	vcftools --vcf /data/home/students/b.huo/new_dataset/filtering_steps/snps_chr.vcf --keep Indica2.txt --max-missing 0.8 --maf 0.25 --recode --out indica2
	vcftools --vcf /data/home/students/b.huo/new_dataset/filtering_steps/snps_chr.vcf --keep Indica.txt --max-missing 0.8 --maf 0.11 --recode --out indica


## 1.2 SMP
#### sort chromosome 
	bcftools sort new_dataset/filtering_steps/meth_chr.vcf  -o new_dataset/filtering_steps/meth_chr.vcf
####  extract indica
	vcftools --vcf /data/home/students/b.huo/new_dataset/filtering_steps/meth_chr.vcf --keep Indica1.txt --max-missing 0.6 --maf 0.2 --recode --out m_indica1
	vcftools --vcf /data/home/students/b.huo/new_dataset/filtering_steps/meth_chr.vcf --keep Indica2.txt --max-missing 0.6 --maf 0.25 --recode --out m_indica2
	vcftools --vcf /data/home/students/b.huo/new_dataset/filtering_steps/meth_chr.vcf --keep Indica.txt --max-missing 0.6 --maf 0.11 --recode --out m_indica

# 2. per site Pi by vcftools
## 2.1 SNP
	cd /data/home/students/b.huo/new_dataset/pops
	vcftools --vcf indica1.recode.vcf --site-pi --out ../pi/indica1
	vcftools --vcf indica2.recode.vcf --site-pi --out ../pi/indica2
	vcftools --vcf indica.recode.vcf --site-pi --out ../pi/indica
#### R
	setwd("~/Desktop/Analysis/new_dataset/calculations/pi")
	Indica1<-read.table("indica1.sites.pi",header=T)
	Indica2<-read.table("indica2.sites.pi",header=T)
	Indica<-read.table("indica.sites.pi",header=T)
	Indica1$indiv<-"Indica I"
	Indica2$indiv<-"Indica II"
	Indica$indiv<-"Indica"
	df<-rbind(Indica1,Indica2,Indica)
	library(ggplot2)
	options(scipen=200)
	pi_snp <- ggplot(df,aes(x=indiv,y=Pi,fill=indiv))+
			  labs(title="Pi_SNP (vcftools)",x="")+
			  geom_boxplot()+
			  ylim(c(0.2,0.9))+
			  guides(fill = guide_legend(title = ""))+
			  scale_fill_manual(values = c("Indica" = "white", "Indica I" = "#F8766D", "Indica II" = "#00BFC4"))+
			  theme_bw()

## 2.2 SMP
	vcftools --vcf m_indica1.recode.vcf --site-pi --out ../pi/m_indica1
	vcftools --vcf m_indica2.recode.vcf --site-pi --out ../pi/m_indica2
	vcftools --vcf m_indica.recode.vcf --site-pi --out ../pi/m_indica

#### R
	setwd("~/Desktop/Analysis/new_dataset/calculations/pi")
	Indica1<-read.table("m_indica1.sites.pi",header=T)
	Indica2<-read.table("m_indica2.sites.pi",header=T)
	Indica<-read.table("m_indica.sites.pi",header=T)
	Indica1$indiv<-"Indica I"
	Indica2$indiv<-"Indica II"
	Indica$indiv<-"Indica"
	df<-rbind(Indica1,Indica2,Indica)
	library(ggplot2)
	options(scipen=200)

	pi_smp <- ggplot(df,aes(x=indiv,y=Pi,fill=indiv))+
			  geom_boxplot()+
			  labs(title="Pi_SMP (vcftools)",x="")+
			  ylim(c(0.2,0.9))+
			  guides(fill = guide_legend(title = ""))+
			  scale_fill_manual(values = c("Indica" = "white", "Indica I" = "#F8766D", "Indica II" = "#00BFC4"))+
			  theme_bw()

	combined <-grid.arrange(pi_snp, pi_smp, ncol = 2)

# 3. PopLDdecay_local machine
	https://github.com/hewm2008/PopLDdecay
	tar -zxvf PopLDdecay-3.42.tar.gz
	cd PopLDdecay-3.42
	cd src
	make
	make clean
	../bin/PopLDdecay
	cd /Users/binghuo/Downloads/PopLDdecay-3.42
## 3.1 SNP	
	# indica1 and indica2 
	./bin/PopLDdecay -InVCF  indica1.recode.vcf.gz -MAF 0.2 -Miss 0.8 -Het 1 -OutStat ./D_300/indica1
	./bin/PopLDdecay -InVCF  indica2.recode.vcf.gz -MAF 0.25 -Miss 0.8 -Het 1 -OutStat ./D_300/indica2
	perl bin/Plot_MultiPop.pl  -inList  ./indica_12.list -keepR -output ./D_300/indica_12
	perl bin/Plot_MultiPop.pl  -inList  ./indica_12.list  -maxX 5 -keepR -output ./plot_D5/indica_12	
	
	# indica 
	./bin/PopLDdecay -InVCF  indica.recode.vcf.gz  -MAF 0.11 -Miss 0.6 -Het 1 -OutStat ./D_300/indica
	perl bin/Plot_OnePop.pl -inFile ./D_300/indica.stat.gz  -keepR -output ~/Downloads/PopLDdecay-3.42/D_300/indica
	perl bin/Plot_OnePop.pl -inFile ./D_300/indica.stat.gz  -maxX 5 -keepR -output ~/Downloads/PopLDdecay-3.42/plot_D5/indica
##### customize in R
## 3.2 SMP
	# indica1 and indica2
	./bin/PopLDdecay -InVCF  m_indica1.recode.vcf.gz -MAF 0.2 -Miss 0.6 -Het 1 -OutStat ./D_300/m_indica1
	./bin/PopLDdecay -InVCF  m_indica2.recode.vcf.gz -MAF 0.25 -Miss 0.6 -Het 1 -OutStat ./D_300/m_indica2
	perl bin/Plot_MultiPop.pl  -inList  ./m_indica_12.list  -keepR -output ~/Downloads/PopLDdecay-3.42/D_300/m_indica_12
	perl bin/Plot_MultiPop.pl  -inList  ./m_indica_12.list -maxX 5 -keepR -output ~/Downloads/PopLDdecay-3.42/plot_D5/m_indica_12
	
	# indica
	./bin/PopLDdecay -InVCF  m_indica.recode.vcf.gz -MAF 0.11 -Miss 0.6 -Het 1 -OutStat ./D_300/m_indica
	perl bin/Plot_OnePop.pl -inFile ./D_300/m_indica.stat.gz  -keepR -output ~/Downloads/PopLDdecay-3.42/D_300/m_indica
	perl bin/Plot_OnePop.pl -inFile ./D_300/m_indica.stat.gz  -maxX 5 -keepR -output ~/Downloads/PopLDdecay-3.42/plot_D5/m_indica

##### find the scripts obtained from keepR and customize in R

# 4. Tajima's D and Pi per gene (Indica)
## 4.1 SNP
#### filtered indels 
	cd new_dataset/snp_consensus
	grep -v "^#" indica.vcf | grep -v "*" > indica1.vcf
	grep -v "^#" indica1.vcf | grep -Ev '.\/\.' > indica.vcf
	grep -v '^#' indica.vcf | wc -l
	56674
	bgzip indica.vcf
	tabix -p vcf indica.vcf.gz

#### create 12 folders named as 1 to 12 and generate consensus
#####	 `-r` is an option that tells `read` to interpret backslashes as literal characters; `-p` option is used to specify an output prefix for the consensus sequence files
	cd /data/home/students/b.huo/new_dataset/snp_consensus
	while read -r value1 value2; do
    while read -r sample; do
        chromosome="$value1" 
        echo "Processing sample: $sample"
        samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$value1:$value2" | bcftools consensus -M N -s "$sample" -p "${sample}_$chromosome_" -H 1pIu indica.vcf.gz >> "$chromosome/${value1}_${value2}_rice.fa"
    done < sample_list
	done < chr_pos
#### R pegas: for loop to calculate tajimaD and Pi for each gene 
	library(ape)
	library(adegenet)
	require(ape)
	library(pegas)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	setwd("/Users/binghuo/Downloads/snp_consensus")
	folder_path="/Users/binghuo/Downloads/snp_consensus"
	snp_indica <- data.frame(Sequence = character(), D = numeric(), Pval.normal = numeric(), Pval.beta = numeric(),
                      Pi=numeric(),Pi_variance=numeric())
	fasta_files <- list.files(folder_path, pattern = "*.fa", full.names = F,recursive=T)

	for (file in fasta_files){
	  x <-read.dna(file,format="fasta",as.character = T)
	  x<-as.DNAbin(x)
	  D_result <- tajima.test(x)
	  Pi_result <-nuc.div(x,TRUE)
	  result_df <- data.frame(Sequence = file,D = D_result$D,Pval.normal = D_result$Pval.normal,Pval.beta = D_result$Pval.beta,Pi=Pi_result[1], Pi_variance=Pi_result[2])
	  next_row <- nrow(snp_indica)+1
	  snp_indica[next_row,] <- result_df
	}
	write.csv(snp_indica,file="SNP_D&pi.csv",quote=FALSE,row.names=FALSE)
## 4.2 SMP 
#### filtered indels 
	cd /data/home/students/b.huo/new_dataset/smp_consensus
	bcftools view -S samples_list m_indica.vcf > m_indica.vcf
	grep -v "^#" m_indica.vcf | grep -v "./\." > m_indica.vcf
	grep -v '^#' m_indica.vcf | wc -l
	1028
	bgzip m_indica.vcf
	tabix -p vcf m_indica.vcf.gz 
#### update REF from reference genome
	bcftools +fill-from-fasta m_indica.vcf -- -c REF -f Oryza_sativa.IRGSP-1.0.dna.toplevel.fa > ref_m_indica.vcf
#### save smp variants sites
	grep -v "^#" m_indica.vcf | cut -f 1,2 > smp_sites
#### extract each base on smp sites from reference sequence
	while read -r chr pos; do
        header_base=$(samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$chr:$pos-$pos")
        base=$(echo "$header_base" | sed -n "2p")
       echo -e "$chr\t$pos\t$base" >> smp_bases   
    done < smp_sites
#### R: change ALT in VCF: most sites of REF are C,G or A, only  6 sites  are "T"--change the corresponding sites with "T" in ALT field as "A" 
	setwd("/Users/binghuo/Downloads/smp_consensus")
	library(dplyr)
	vcf_file <- "ref_m_indica.vcf"
	smp_bases_file <- "smp_bases"
	smp_bases_df <- read.table(smp_bases_file, header = FALSE, col.names = c("chr", "pos", "base"), sep = "\t")
	
	###### Filter smp_bases for base "T"
	smp_bases_filtered <- smp_bases_df %>% filter(base == "T")
	
	###### Create a lookup table for positions to replace with "A"
	positions_to_replace <- as.data.frame(cbind(smp_bases_filtered$chr, smp_bases_filtered$pos))
	colnames(positions_to_replace)<-c("chr","pos")
	
	###### Read the VCF file and replace
	vcf_df <- read.table(vcf_file, header = FALSE, col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT","C148","C019","W169","MH63","C151","C139","C135","W161","ZS97"), sep = "\t",stringsAsFactors = FALSE, colClasses = c("character"))                                                 
	vcf_df$ALT[vcf_df$CHROM %in% positions_to_replace$chr & vcf_df$POS %in% positions_to_replace$pos] <- "A"
	write.table(vcf_df, file = "final_m_indica.vcf", quote = FALSE, sep = "\t", row.names = FALSE)

	tabix -p vcf final_m_indica.vcf.gz
#### create 12 folders named as 1 to 12 and generate consensus
	while read -r value1 value2; do
    while read -r sample; do
        chromosome="$value1" 
        echo "Processing sample: $sample"
        samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$value1:$value2" | bcftools consensus -M N -s "$sample" -p "${sample}_$chromosome_" -H 1pIu final_m_indica.vcf.gz >> "$chromosome/${value1}_${value2}_rice.fa"
    done < sample_list
	done < chr_pos
#### R pegas: for loop to calculate tajimaD and Pi for each gene 
	library(ape)
	library(adegenet)
	require(ape)
	library(pegas)
	library(dplyr)
	library(ggplot2)
	library(gridExtra)
	setwd("/Users/binghuo/Downloads/smp_consensus")
	folder_path="/Users/binghuo/Downloads/smp_consensus"
	smp_indica <- data.frame(Sequence = character(), D = numeric(), Pval.normal = numeric(), Pval.beta = numeric(),Pi=numeric(),Pi_variance=numeric())
	fasta_files <- list.files(folder_path, pattern = "*.fa", full.names = F,recursive=T)


	for (file in fasta_files){
	  y <-read.dna(file,format="fasta",as.character = T)
	  y <-as.DNAbin(y)
	  D_result <- tajima.test(y)
	  Pi_result <-nuc.div(y,TRUE)
	  result_df <- data.frame(Sequence = file,D = D_result$D,Pval.normal = D_result$Pval.normal,Pval.beta = D_result$Pval.beta,
                          Pi=Pi_result[1], Pi_variance=Pi_result[2])
	  next_row <- nrow(smp_indica)+1
	  smp_indica[next_row,] <- result_df
	}
	write.csv(smp_indica,file="SMP_D&pi.csv",quote=FALSE,row.names=FALSE)

	D1<-ggplot()+
	  geom_boxplot(mapping=aes(x="Indica",y=snp_indica$D))+
	  labs(title="Tajima's D (SNP)",x="",y="Tajima.D")+
	  ylim(c(-4,3))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()
	D2<-ggplot()+
	  geom_boxplot(mapping=aes(x="Indica",y=smp_indica$D))+
	  labs(title="Tajima's D (SMP)",x="",y="Tajima.D")+
	  ylim(c(-4,3))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()

	f_smp_indica <-smp_indica %>% filter(Pi >0)
	f_snp_indica <-snp_indica %>% filter(Pi >0)
	min_value <- min(min(f_snp_indica$Pi),min(f_smp_indica$Pi))
	max_value <- max(max(f_snp_indica$Pi), 	max(f_smp_indica$Pi))

	P1<- ggplot()+
	  geom_boxplot(mapping=aes(x="Indica",y=snp_indica$Pi))+
	  labs(title="Pi_log10 scale (SNP)",x="",y="Pi_log10 scale")+
	  scale_y_log10(limits=c(min_value,max_value))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()
	P2<-ggplot()+
	  geom_boxplot(mapping=aes(x="Indica",y=smp_indica$Pi))+
	  labs(title="Pi_log10 scale (SMP)",x="",y="Pi_log10 scale")+
	  scale_y_log10(limits=c(min_value,max_value))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()

	combined_D <-grid.arrange(D1, D2, ncol = 2)
	combined_Pi <-grid.arrange(P1, P2, ncol = 2)
# 5. Tajima's D and Pi per gene (Indica I and Indica II)
## 5.1 SNP
#### filtered indels 
	cd new_dataset/snp_consensus
	grep -v "^#" indica.vcf | grep -v "*" > indica1.vcf
	grep -v "^#" indica1.vcf | grep -Ev '.\/\.' > indica.vcf
	grep -v '^#' indica.vcf | wc -l
	56674
	bgzip indica.vcf
	tabix -p vcf indica.vcf.gz

#### create 12 folders named as 1 to 12 and generate consensus
#####	 `-r` is an option that tells `read` to interpret backslashes as literal characters; `-p` option is used to specify an output prefix for the consensus sequence files
	cd /data/home/students/b.huo/new_dataset/snp_consensus/snp_indica1
	bcftools view -s "C019,C135,C139,C151,ZS97" -o indica1.vcf indica.vcf
	bgzip indica1.vcf
	tabix -p vcf indica1.vcf.gz
	while read -r value1 value2; do
    while read -r sample; do
        chromosome="$value1" 
        echo "Processing sample: $sample"
        samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$value1:$value2" | bcftools consensus -M N -s "$sample" -p "${sample}_$chromosome_" -H 1pIu indica1.vcf.gz >> "$chromosome/${value1}_${value2}_rice.fa"
    done < indica1.txt
	done < chr_pos
	
	cd /data/home/students/b.huo/new_dataset/snp_consensus/snp_indica2
	bcftools view -s "C148,W169,MH63,W161" -o indica2.vcf indica.vcf
	bgzip indica2.vcf
	tabix -p vcf indica2.vcf.gz
	while read -r value1 value2; do
    while read -r sample; do
        chromosome="$value1" 
        echo "Processing sample: $sample"
        samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$value1:$value2" | bcftools consensus -M N -s "$sample" -p "${sample}_$chromosome_" -H 1pIu indica2.vcf.gz >> "$chromosome/${value1}_${value2}_rice.fa"
    done < indica2.txt
	done < chr_pos
#### R pegas: for loop to calculate tajimaD and Pi for each gene 
	library(ape)
	library(adegenet)
	require(ape)
	library(pegas)
	library(ggplot2)
	library(dplyr)
	library(gridExtra)
	setwd("/Users/binghuo/Downloads/snp_consensus/snp_indica1")
	folder_path="/Users/binghuo/Downloads/snp_consensus/snp_indica1"
	snp_indica1 <- data.frame(Sequence = character(), D = numeric(), Pval.normal = numeric(), Pval.beta = numeric(),Pi=numeric(),Pi_variance=numeric())
	fasta_files <- list.files(folder_path, pattern = "*.fa", full.names = F,recursive=T)

	for (file in fasta_files){
	  x <-read.dna(file,format="fasta",as.character = T)
	  x<-as.DNAbin(x)
	  D_result <- tajima.test(x)
	  Pi_result <-nuc.div(x,TRUE)
	  result_df <- data.frame(Sequence = file,D = D_result$D,Pval.normal = D_result$Pval.normal,Pval.beta = D_result$Pval.beta,
                          Pi=Pi_result[1], Pi_variance=Pi_result[2])
	  next_row <- nrow(snp_indica1)+1
	  snp_indica1[next_row,] <- result_df
	}
	write.csv(snp_indica1,file="SNP_indica1_D&pi.csv",quote=FALSE,row.names=FALSE)

	#######################
	setwd("/Users/binghuo/Downloads/snp_consensus/snp_indica2")
	folder_path="/Users/binghuo/Downloads/snp_consensus/snp_indica2"
	snp_indica2 <- data.frame(Sequence = character(), D = numeric(), Pval.normal = numeric(), Pval.beta = numeric(),Pi=numeric(),Pi_variance=numeric())
	fasta_files <- list.files(folder_path, pattern = "*.fa", full.names = F,recursive=T)
	for (file in fasta_files){
	  x <-read.dna(file,format="fasta",as.character = T)
	  x<-as.DNAbin(x)
	  D_result <- tajima.test(x)
	  Pi_result <-nuc.div(x,TRUE)
	  result_df <- data.frame(Sequence = file,D = D_result$D,Pval.normal = D_result$Pval.normal,Pval.beta = D_result$Pval.beta,
                          Pi=Pi_result[1], Pi_variance=Pi_result[2])
	  next_row <- nrow(snp_indica2)+1
	  snp_indica2[next_row,] <- result_df
	}
	write.csv(snp_indica2,file="SNP_indica2_D&pi.csv",quote=FALSE,row.names=FALSE)

## 5.2 SMP 
#### filtered indels 
	cd /data/home/students/b.huo/new_dataset/smp_consensus
	bcftools view -S samples_list m_indica.vcf > m_indica.vcf
	grep -v "^#" m_indica.vcf | grep -v "./\." > m_indica.vcf
	grep -v '^#' m_indica.vcf | wc -l
	1028
	bgzip m_indica.vcf
	tabix -p vcf m_indica.vcf.gz 
#### update REF from reference genome
	bcftools +fill-from-fasta m_indica.vcf -- -c REF -f Oryza_sativa.IRGSP-1.0.dna.toplevel.fa > ref_m_indica.vcf
#### save smp variants sites
	grep -v "^#" m_indica.vcf | cut -f 1,2 > smp_sites
#### extract each base on smp sites from reference sequence
	while read -r chr pos; do
        header_base=$(samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$chr:$pos-$pos")
        base=$(echo "$header_base" | sed -n "2p")
       echo -e "$chr\t$pos\t$base" >> smp_bases   
    done < smp_sites
#### R: change ALT in VCF: most sites of REF are C,G or A, only  6 sites  are "T"--change the corresponding sites with "T" in ALT field as "A" 
	setwd("/Users/binghuo/Downloads/smp_consensus")
	library(dplyr)
	vcf_file <- "ref_m_indica.vcf"
	smp_bases_file <- "smp_bases"
	smp_bases_df <- read.table(smp_bases_file, header = FALSE, col.names = c("chr", "pos", "base"), sep = "\t")
	
	###### Filter smp_bases for base "T"
	smp_bases_filtered <- smp_bases_df %>% filter(base == "T")
	
	###### Create a lookup table for positions to replace with "A"
	positions_to_replace <- as.data.frame(cbind(smp_bases_filtered$chr, smp_bases_filtered$pos))
	colnames(positions_to_replace)<-c("chr","pos")
	
	###### Read the VCF file and replace
	vcf_df <- read.table(vcf_file, header = FALSE, col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT","C148","C019","W169","MH63","C151","C139","C135","W161","ZS97"), sep = "\t",stringsAsFactors = FALSE, colClasses = c("character"))                                                 
	vcf_df$ALT[vcf_df$CHROM %in% positions_to_replace$chr & vcf_df$POS %in% positions_to_replace$pos] <- "A"
	write.table(vcf_df, file = "final_m_indica.vcf", quote = FALSE, sep = "\t", row.names = FALSE)

	tabix -p vcf final_m_indica.vcf.gz
#### create 12 folders named as 1 to 12 and generate consensus
	cd /data/home/students/b.huo/new_dataset/smp_consensus/smp_indica1
	bcftools view -s "C019,C135,C139,C151,ZS97" -o m_indica1.vcf final_m_indica.vcf
	bgzip m_indica1.vcf
	tabix -p vcf m_indica1.vcf.gz
	while read -r value1 value2; do
    while read -r sample; do
        chromosome="$value1" 
        echo "Processing sample: $sample"
        samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$value1:$value2" | bcftools consensus -M N -s "$sample" -p "${sample}_$chromosome_" -H 1pIu m_indica1.vcf.gz >> "$chromosome/${value1}_${value2}_rice.fa"
    done < indica1.txt
	done < chr_pos
	
	cd /data/home/students/b.huo/new_dataset/smp_consensus/smp_indica2
	bcftools view -s "C148,W169,MH63,W161" -o m_indica2.vcf final_m_indica.vcf
	bgzip m_indica2.vcf
	tabix -p vcf m_indica2.vcf.gz
	while read -r value1 value2; do
    while read -r sample; do
        chromosome="$value1" 
        echo "Processing sample: $sample"
        samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa "$value1:$value2" | bcftools consensus -M N -s "$sample" -p "${sample}_$chromosome_" -H 1pIu m_indica2.vcf.gz >> "$chromosome/${value1}_${value2}_rice.fa"
    done < indica2.txt
	done < chr_pos
#### R pegas: for loop to calculate tajimaD and Pi for each gene 
	library(ape)
	library(adegenet)
	require(ape)
	library(pegas)
	library(ggplot2)
	library(dplyr)
	library(gridExtra)
	setwd("/Users/binghuo/Downloads/smp_consensus/smp_indica1")
	folder_path="/Users/binghuo/Downloads/smp_consensus/smp_indica1"
	smp_indica1 <- data.frame(Sequence = character(), D = numeric(), Pval.normal = numeric(), Pval.beta = numeric(),Pi=numeric(),Pi_variance=numeric())
	fasta_files <- list.files(folder_path, pattern = "*.fa", full.names = F,recursive=T)

	for (file in fasta_files){
	  x <-read.dna(file,format="fasta",as.character = T)
	  x<-as.DNAbin(x)
	  D_result <- tajima.test(x)
	  Pi_result <-nuc.div(x,TRUE)
	  result_df <- data.frame(Sequence = file,D = D_result$D,Pval.normal = D_result$Pval.normal,Pval.beta = D_result$Pval.beta,
                          Pi=Pi_result[1], Pi_variance=Pi_result[2])
	  next_row <- nrow(smp_indica1)+1
	  smp_indica1[next_row,] <- result_df
	}
	write.csv(smp_indica1,file="SMP_indica1_D&pi.csv",quote=FALSE,row.names=FALSE)

	#######################
	setwd("/Users/binghuo/Downloads/smp_consensus/smp_indica2")
	folder_path="/Users/binghuo/Downloads/smp_consensus/smp_indica2"
	smp_indica2 <- data.frame(Sequence = character(), D = numeric(), Pval.normal = numeric(), Pval.beta = numeric(),Pi=numeric(),Pi_variance=numeric())
	fasta_files <- list.files(folder_path, pattern = "*.fa", full.names = F,recursive=T)
	for (file in fasta_files){
	  x <-read.dna(file,format="fasta",as.character = T)
	  x<-as.DNAbin(x)
	  D_result <- tajima.test(x)
	  Pi_result <-nuc.div(x,TRUE)
	  result_df <- data.frame(Sequence = file,D = D_result$D,Pval.normal = D_result$Pval.normal,Pval.beta = D_result$Pval.beta,
                          Pi=Pi_result[1], Pi_variance=Pi_result[2])
	  next_row <- nrow(smp_indica2)+1
	  smp_indica2[next_row,] <- result_df
	}
	write.csv(smp_indica2,file="SMP_indica2_D&pi.csv",quote=FALSE,row.names=FALSE)

	d1<- ggplot()+
	  geom_boxplot(mapping=aes(x="Indica I",y=snp_indica1$D,fill="Indica I"))+
	  geom_boxplot(mapping=aes(x="Indica II",y=snp_indica2$D,fill="Indica II"))+
	  labs(title="Tajima's D (SNP)",x="",y="Tajima.D")+
	  ylim(c(-10,5))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()

	d2 <- ggplot()+
	  geom_boxplot(mapping=aes(x="Indica I",y=smp_indica1$D,fill="Indica I"))+
	  geom_boxplot(mapping=aes(x="Indica II",y=smp_indica2$D,fill="Indica II"))+
	  labs(title="Tajima's D (SMP)",x="",y="Tajima.D")+
	  ylim(c(-10,5))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()

	f_snp_indica1 <-snp_indica1 %>% filter(Pi >0)
	f_snp_indica2 <-snp_indica2 %>% filter(Pi >0)
	f_smp_indica1 <-smp_indica1 %>% filter(Pi >0)
	f_smp_indica2 <-smp_indica2 %>% filter(Pi >0)
	min_value <- min(min(f_snp_indica1$Pi), min(f_snp_indica2$Pi),min(f_smp_indica1$Pi),min(f_smp_indica2$Pi))
	max_value <- max(max(f_snp_indica1$Pi), max(f_snp_indica2$Pi),max(f_smp_indica1$Pi),max(f_smp_indica2$Pi))

	p1<- ggplot()+
	  geom_boxplot(mapping=aes(x="Indica I",y=f_snp_indica1$Pi,fill="Indica I"))+
	  geom_boxplot(mapping=aes(x="Indica II",y=f_snp_indica2$Pi,fill="Indica II"))+
	  labs(title="Pi_log10 scale (SNP)",x="",y="Pi_log10 scale")+
	  scale_y_log10(limits=c(min_value,max_value))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()
	p2 <- ggplot()+
	  geom_boxplot(mapping=aes(x="Indica I",y=f_smp_indica1$Pi,fill="Indica I"))+
	  geom_boxplot(mapping=aes(x="Indica II",y=f_smp_indica2$Pi,fill="Indica II"))+
	  labs(title="Pi_log10 scale (SMP)",x="",y="Pi_log10 scale")+
	  scale_y_log10(limits=c(min_value,max_value))+
	  guides(fill = guide_legend(title = ""))+
	  theme_bw()

	combined_d <-grid.arrange(d1, d2, ncol = 2)
	combined_pi <-grid.arrange(p1, p2, ncol = 2)
