# Filtering steps

### SNP
	grep -v '^#' rice_snps_filtered.recode.vcf | wc -l
##### *rice_snps_filtered.recode.vcf:     1,251,872 *

### 1. remove non_numerical chromosomes and keep biallelic sites
	sed '/#/d' ./rice_snps_filtered.recode.vcf | cut -f 1 | sort -u
	sed -e '/^[A-Z]/d' rice_snps_filtered.recode.vcf> ./new_dataset/snps_chr.vcf
	sed '/#/d' ./new_dataset/snps_chr.vcf | cut -f 1 | sort -u
	
	grep -v '^#' ./new_dataset/snps_chr.vcf | wc -l
##### * snps_chr.vcf:     1, 247,596*
###### **- m2: only view biallelic SNPs**
	bcftools view -v snps -m2 ./new_dataset/snps_chr.vcf > ./new_dataset/snps_biallelic.vcf
	grep -v  '^#' ./new_dataset/snps_biallelic.vcf | wc -l
##### * snps_biallelic.vcf:     1,247,596*

### 2. rm intermediate 
	bcftools view -s "C019,C135,C139,C151,ZS97,C148,W161,W169,MH63,W105,W286,W081,W306,W257,NIP" -o ./new_dataset/snps_chr_rm.vcf ./new_dataset/snps_biallelic.vcf	
	grep -v  '^#'  ./new_dataset/snps_chr_rm.vcf | wc -l
##### *snps_chr_rm.vcf :     1,247,596 ( 82.94% NA)

### 3. filtering NA 80%

	vcftools --vcf ./new_dataset/snps_chr_rm.vcf  --max-missing 0.8 --recode --out ./new_dataset/snps_chr_rm_na
	grep -v  '^#' ./new_dataset/snps_chr_rm_na.recode.vcf | wc -l
##### *snps_chr_rm_na.recode.vcf:     144, 787 ( 7.93% NA)
###### LD 100 5 0.2,    --geno 0.1    
	plink --vcf ./new_dataset/snps_chr_rm_na.recode.vcf  --indep-pairwise 100 5 0.2 --geno 0.1 --maf 0.05 --hwe 0.0001 --make-bed --out ./new_dataset/new_admixture/snp --allow-extra-chr --set-missing-var-ids @:#
	grep -v '^#' ./new_dataset/new_admixture/snp.bim | wc -l

###### *snp.bim 80,520
	plink --bfile ./new_dataset/new_admixture/snp --out ./new_dataset/new_admixture/snp_for_vcf --recode vcf
	plink --vcf ./new_dataset/new_admixture/snp_for_vcf.vcf --extract ./new_dataset/new_admixture/snp.prune.out --out ./new_dataset/new_admixture/snp_filtered --recode vcf-iid --keep-allele-order --allow-extra-chr --set-missing-var-ids @:#
	grep -v '^#' ./new_dataset/new_admixture/snp_filtered.vcf | wc -l
#####  * snp_filtered.vcf  78,341 

	plink --vcf ./new_dataset/new_admixture/snp_filtered.vcf  --make-bed --out ./new_dataset/new_admixture/snp_filtered
##### *snp_filtered.bed  78,341

### Methylation

	grep -v '^#' rice_meth_all.vcf | wc -l
##### * rice_meth_all.vcf:      28,411*

### 1. rm non_numerical chromosomes 
	sed '/#/d' ./rice_meth_all.vcf | cut -f 1 | sort -u
	sed -e '/^[A-Z]/d' rice_meth_all.vcf> ./new_dataset/meth_chr.vcf
	sed '/#/d' ./new_dataset/meth_chr.vcf | cut -f 1 | sort -u
	grep -v '^#' ./new_dataset/meth_chr.vcf | wc -l
##### * meth_chr.vcf:     18,197 *
	
### 2. rm intermediate 
	bcftools view -s "C019,C135,C139,C151,ZS97,C148,W161,W169,MH63,W105,W286,W081,W306,W257,NIP" -o ./new_dataset/meth_chr_rm.vcf ./new_dataset/meth_chr.vcf
	grep -v  '^#'  ./new_dataset/meth_chr_rm.vcf | wc -l
##### *meth_chr_rm.vcf :     18,197 * (34.53 % NA)
#####  sort	the positions in order
	bcftools sort ./new_dataset/meth_chr_rm.vcf  -o output_sorted.vcf
##### change *output_sorted.vcf* to *meth_chr_rm.vcf*

### 3. filtering NA 80%

	vcftools --vcf ./new_dataset/meth_chr_rm.vcf  --max-missing 0.8 --recode --out ./new_dataset/meth_chr_rm_na
	grep -v  '^#' ./new_dataset/meth_chr_rm_na.recode.vcf | wc -l
##### *meth_chr_rm_na.recode.vcf* :   3,766  (9.81% NA)
	plink --vcf ./new_dataset/meth_chr_rm_na.recode.vcf --indep-pairwise 100 5 0.2 --maf 0.05 --hwe 0.0001 --make-bed --out ./new_dataset/new_admixture/meth --allow-extra-chr --set-missing-var-ids @:#
	plink --bfile ./new_dataset/new_admixture/meth --out ./new_dataset/meth_for_vcf --recode vcf
	grep -v '^#' ./new_dataset/new_admixture/meth.bim  | wc -l 
##### *meth.bim*: 3,253 
	plink --vcf ./new_dataset/meth_for_vcf.vcf --extract ./new_dataset/new_admixture/meth.prune.out --out ./new_dataset/new_admixture/meth_filtered  --recode vcf-iid --keep-allele-order --allow-extra-chr --set-missing-var-ids @:#
	grep -v '^#' ./new_dataset/new_admixture/meth_filtered.vcf  | wc -l 
##### * meth_filtered.vcf*:    3,034 
	plink --vcf ./new_dataset/new_admixture/meth_filtered.vcf   --make-bed --out ./new_dataset/new_admixture/meth_filtered
##### * meth_filtered.bed* :   3,034



# Phylogeny tree

### SNP 
##### *snp_filtered.vcf*:  78,341 
##### downloads VCF2Dis first
	cd downloads
	cd VCF2Dis-1.50
	 ./bin/VCF2Dis -InPut  ./bin/snp_filtered.vcf  -OutPut snp_p_dis.mat
	 
http://www.atgc-montpellier.fr/fastme/
##### Choose input data:  mat file, distance matrix; tree building: NJ// download nwk file

####  R tree visualization and customize
	library(ggtree)
	library(treeio)
	library(ggplot2)
	library(dplyr)
	library(ape)
	setwd("~/Desktop/Analysis/new_dataset/phylogeny_new")
	snp_tree<-read.tree("snp_p_dis_mat_fastme-tree.nwk")
	tip.group <- read.csv('sample_group_phylogeny.txt',sep=";",header=TRUE)
	tip.group
	tree.a1 <- full_join(snp_tree,tip.group,by="label")
	ggtree(tree.a1, 
       layout = "ape",
       size=0.5,
       aes(color=group))+
	geom_tiplab(size=4,face="italic")+
	  xlim(-0.1,0.5)+
	  ylim(-0.1,0.5)+
	  labs(title ="SNP")+
	  theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.55, vjust = 0.1),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))+
	  guides(color = guide_legend(override.aes = list(label = "")))+
	  geom_treescale(x = -0.1, y = -0.1, offset = 0.02, fontsize = 3)+
	  scale_color_manual(values = c("Indica I" = "#00BA38", "Indica II" = "#00BFC4", "Aus" = "#B79F00",
                                "Aromatic"="#F8766D","Japonica"="#F564E3"))

### Methylation
##### *meth_filtered.bed*:   3,034 
	cd downloads
	cd VCF2Dis-1.50
	 ./bin/VCF2Dis -InPut  ./bin/meth_filtered.vcf  -OutPut meth_p_dis.mat

http://www.atgc-montpellier.fr/fastme/
#### Choose input data:  mat file, distance matrix; tree building: NJ// download nwk file

####  R tree visualization and customize
	library(ggtree)
	library(treeio)
	library(ggplot2)
	library(dplyr)
	library(ape)
	setwd("~/Desktop/Analysis/new_dataset/phylogeny_new")
	meth_tree<-read.tree("meth_p_dis_mat_fastme-tree.nwk")
	tip.group <- read.csv('sample_group_phylogeny.txt',sep=";",header=TRUE)
	tip.group
	tree.a1 <- full_join(meth_tree,tip.group,by="label")
	ggtree(tree.a1, 
       layout = "ape",
       size=0.5,
       aes(color=group))+
	  geom_tiplab(size=4,face="italic")+
	  xlim(-0.1,0.4)+
	  ylim(-0.1,0.4)+
	  labs(title ="Meth")+
	  theme(legend.text = element_text(face="italic",size=10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.4, vjust = -4),
        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"))+
	  guides(color = guide_legend(override.aes = list(label = "")))+
	  geom_treescale(x = -0.1, y = -0.1, offset = 0.02, fontsize = 3)+
	  scale_color_manual(values = c("Indica I" = "#00BA38", "Indica II" = "#00BFC4", "Aus" = "#B79F00",
                                "Aromatic"="#F8766D","Japonica"="#F564E3"))
	

# Admxiture
### SNP  
##### *snp_filtered.bed*:   78,341
	cd ./new_dataset/new_admixture
	for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ; do admixture --cv snp_filtered.bed $K|tee snp_filtered_log${K}.out; done
	grep -h CV snp_filtered_log*.out
##### save CV
	grep -h CV snp_filtered_log*.out|awk -F ':' '{print NR"\t"$2}'|sed '1i\\K\tCV_error' > snp_CV_for_plot.txt 
#### plot CV values in R , choose k for structural plots
	library(ggplot2)
	setwd("~/Desktop/Analysis/new_dataset/admixture")
	mydata<-read.table("snp_CV_for_plot.txt",header=T,sep="\t")
	qplot(x=K,y=CV_error,color=I('black'),data=mydata)+geom_line(color="red",lwd=1)+
	  ggtitle('snp_CV')+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))	
##### upload the Q files and plot pophelper.com 

### Methylation
##### *meth_filtered.bed*:    3,034
	cd ./new_dataset/new_admixture
	for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ; do admixture --cv meth_filtered.bed  $K | tee meth_log${K}.out; done
	grep -h CV meth_log*.out
	grep -h CV meth_log*.out|awk -F ':' '{print NR"\t"$2}'|sed '1i\\K\tCV_error' > meth_CV_for_plot.txt 

	library(ggplot2)
	setwd("~/Desktop/Analysis/new_dataset/admixture")
	mydata1<-read.table("meth_CV_for_plot.txt",header=T,sep="\t")
	qplot(x=K,y=CV_error,color=I('black'),data=mydata1)+geom_line(color="red",lwd=1)+
	  ggtitle('meth_CV')+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

##### upload the Q files and plot pophelper.com 


# DAPC (mean imputation)

### SNP 
##### *snp_filtered.vcf*:  78,341 
	library(ape)
	library(pegas)
	library(seqinr)
	library(ggplot2)
	library(adegenet)
	library(gtools)
	library(vcfR)
	library(dplyr)
	library(dartR)

	setwd("~/Desktop/Analysis/new_dataset")
	snp_vcf <- read.vcfR("snp_filtered.vcf")
	snp_genlight <- vcfR2genlight(snp_vcf)
	popNames(snp_genlight)

##### Data frame with group names corresponding to each sample
	sample_data <- data.frame(sample_names = c('C019','C135','C139','C151','ZS97',
                                           'C148','W161','W169','MH63',
                                           'W105','W286',
                                           'W081','W306',
                                           'W257','NIP'))
	group_data <- data.frame(sample_names = c('C019','C135','C139','C151','ZS97',
                                          'C148','W161','W169','MH63',
                                          'W105','W286',
                                          'W081','W306',
                                          'W257','NIP'),
                         Group = c("Indica I","Indica I","Indica I","Indica I","Indica I",
                                   "Indica II","Indica II","Indica II","Indica II",
                                   "Aus","Aus",
                                   "Aromatic","Aromatic",
                                   "Japonica","Japonica"))

	joined_data <- full_join(sample_data, group_data, by = "sample_names")
	group_names <- c("Aromatic","Aus","Indica I","Indica II","Japonica")
	group_colors <- c("#FB6496","#99600F","#2BCE48","#005C31","#C814FA")
	myCols <- setNames(group_colors,group_names)

#### define the group for each individual
	gl_1 <- gl.define.pop(snp_genlight, ind.list=c("C019","C135","C139","C151","ZS97"), new="Indica I", verbose = NULL)
	gl_2 <- gl.define.pop(gl_1, ind.list=c("C148","W161","W169","MH63"), new="Indica II", verbose = NULL)
	gl_3 <- gl.define.pop(gl_2, ind.list=c("W105","W286"), new="Aus", verbose = NULL)
	gl_4 <- gl.define.pop(gl_3, ind.list=c("W081","W306"), new="Aromatic", verbose = NULL)
	gl <- gl.define.pop(gl_4, ind.list=c("W257","NIP"), new="Japonica", verbose = NULL)
	popNames(gl)

#### find.clusters--> k=2; cross validation to choose the best number of PCs
	grp <- find.clusters(gl)
	mat <- as.matrix(gl)
	genotype_matrix <- apply(mat,2,na.replace,mean, na.rm = TRUE)
	xval <- xvalDapc(genotype_matrix, gl$pop,n.pca.max=15,training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE)
##### xval$`Number of PCs Achieving Lowest MSE` ---4
	xval[2:6] 
#### DAPC plot and customize
	dapc1 <- dapc(gl, n.pca=4,n.da=4)

	scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=myCols, solid=.4, cex=5, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.leg ="bottomleft",
        leg=TRUE, txt.leg=group_names)
	title("SNP_DAPC plot",font.main=1)
	par(xpd=TRUE,new=TRUE)
	df <- data.frame(x=dapc1$ind.coord[,1],y=dapc1$ind.coord[,2])
	s.label(dfxy=df,xax=1,yax=2,label=group_data$sample_names,
        clabel=0.5,
        grid = FALSE,addaxes = FALSE,boxes=FALSE)

	points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=1, lwd=2, col="black")
	points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=2, lwd=2, col=myCols)

	myInset <- function(){
	temp <- dapc1$pca.eig
	temp <- 100* cumsum(temp)/sum(temp)
	plot(temp, col=rep(c("black","lightgrey"),
                   c(dapc1$n.pca,1000)), ylim=c(0,100),
     xlab="PCA axis", ylab="Cumulated variance (%)",
     cex=1, pch=20, type="h", lwd=2)}
	
	add.scatter(myInset(), posi="topleft",
            inset=c(-0.03,-0.04), ratio=.28,
            bg=transp("white"))

### Methylation 
##### *meth_filtered.bed*:   3,034 
	library(ape)
	library(pegas)
	library(seqinr)
	library(ggplot2)
	library(adegenet)
	library(gtools)
	library(vcfR)
	library(dplyr)
	library(dartR)
	setwd("~/Desktop/Analysis/new_dataset")
	meth_vcf <- read.vcfR("meth_filtered.vcf")
	meth_genlight <- vcfR2genlight(meth_vcf )

	popNames(meth_genlight)

	sample_data <- data.frame(sample_names = c('C019','C135','C139','C151','ZS97',
                                           'C148','W161','W169','MH63',
                                           'W105','W286',
                                           'W081','W306',
                                           'W257','NIP'))

	group_data <- data.frame(sample_names = c('C019','C135','C139','C151','ZS97',
                                          'C148','W161','W169','MH63',
                                          'W105','W286',
                                          'W081','W306',
                                          'W257','NIP'),
                         Group = c("Indica I","Indica I","Indica I","Indica I","Indica I",
                                   "Indica II","Indica II","Indica II","Indica II",
                                   "Aus","Aus",
                                   "Aromatic","Aromatic",
                                   "Japonica","Japonica"))

	joined_data <- full_join(sample_data, group_data, by = "sample_names")
	group_names <- c("Aromatic","Aus","Indica I","Indica II","Japonica")
	group_colors <- c("#FB6496","#99600F","#2BCE48","#005C31","#C814FA")
	myCols <- setNames(group_colors,group_names)
###
	meth_1 <- gl.define.pop(meth_genlight, ind.list=c("C019","C135","C139","C151","ZS97"), new="Indica I", verbose = NULL)
	meth_2 <- gl.define.pop(meth_1, ind.list=c("C148","W161","W169","MH63"), new="Indica II", verbose = NULL)
	meth_3 <- gl.define.pop(meth_2, ind.list=c("W105","W286"), new="Aus", verbose = NULL)
	meth_4 <- gl.define.pop(meth_3, ind.list=c("W081","W306"), new="Aromatic", verbose = NULL)
	meth <- gl.define.pop(meth_4, ind.list=c("W257","NIP"), new="Japonica", verbose = NULL)
###
	grp1 <- find.clusters(meth)
	popNames(meth)
	mat2 <- as.matrix(meth)
	genotype_matrix2 <- apply(mat2,2,na.replace,mean, na.rm = TRUE)
	xval2 <- xvalDapc(genotype_matrix2, meth$pop, n.pca.max = 15, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE)

##### xval2[2:6] #xval2$`Number of PCs Achieving Lowest MSE` --1 

	dapc2 <- dapc(meth,n.pca=1)
	scatter(dapc2, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=myCols, solid=.4, cex=5, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.leg ="topleft",
        leg=TRUE, txt.leg=group_names)
	title("Meth_DAPC plot",font.main=1)
	par(xpd=TRUE,new=TRUE)
	df <- data.frame(x=dapc2$ind.coord[,1],y=dapc2$ind.coord[,2])
	s.label(dfxy=df,xax=1,yax=2,label=group_data$sample_names,
        clabel=0.5,
        grid = FALSE,addaxes = FALSE,boxes=FALSE)
	points(dapc2$grp.coord[,1], 	dapc2$grp.coord[,2], pch=4,
       cex=1, lwd=2, col="black")
	points(dapc2$grp.coord[,1], dapc2$grp.coord[,2], pch=4,
       cex=2, lwd=2, col=myCols)

	myInset <- function(){
	  temp <- dapc2$pca.eig
	  temp <- 100* cumsum(temp)/sum(temp)
	  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc2$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=20, type="h", lwd=2)}
	add.scatter(myInset(), posi="topright",
            inset=c(0,0.04), ratio=0.18,
            bg=transp("white"))

