library(plyr)
library(dplyr)
library(ggplot2)
library(mclust)
library(qtl)

####################################
#remove this before final packaging#
####################################
setwd("~/HTA_Linkage/")
####################################
#remove this before final packaging#
####################################

#Read in the table with Matt Rockmans SNP set
rockmanData <- read.table("GenotypeProcessing/RILS.txt", header=TRUE)
rockmanSNPs <- rockmanData$SNP_ID

# Read in the raw GoldenGate genotyping files

load("GenotypeProcessing/rawplates.RData")

# Put in strain assignments for the raw data from Rockefeller

strains <- read.csv("GenotypeProcessing/GG_strains.csv", header=F, stringsAsFactors=F)[,1]

plate1$strains <- rep(strains[1:96], each=1536)
plate2$strains <- rep(strains[97:192], each=1536)
plate3$strains <- rep(strains[289:384], each=1536)
plate4$strains <- rep(strains[193:288], each=1536)

# Note, plate 3 and plate 4 assignments are switched above to fix the plate mixup at Rockefeller. Can only
# really see mixups with plate 4 because there are non-RIAILs on that plate. Plates 1 and 2 are correct.

plate1$plate <- 1
plate2$plate <- 2
plate3$plate <- 3
plate4$plate <- 4

all <- rbind(plate1, plate2, plate3, plate4)

# Reduce to only useful data (strains, SNP names, positions, and raw genotype data)

gg <- with(all, data.frame(strain = strains, snp=SNP.Name, chr=Chr, pos=Position, x=X.Raw, y=Y.Raw))

# Make an error vector for the tryCatch (sets all classification values to NA)

mclust.error <- Mclust(data=gg[gg$snp=="CE1-28", c("x","y")], G=2, modelNames=c("VVV"), warn=T)

mclust.error$classification[1:384] <- rep(NA, 384)

# Assign genotypes to raw data using EM algorithm. This step takes a few minutes.

all.clusters <- dlply(.data=gg, .variables=.(snp), .fun=function(D) {
    tryCatch({out <- Mclust(data=D[,c("x","y")], G=2, modelNames=c("VVV"), warn=T)
              #out$data<-D[,c("strain", "snp", "chr", "pos")]
              return(out)}, error=function(x){print(D$snp[1]); return(mclust.error)})
})

# Turn the values into numeric for the two values of "z" and that of classification

z1 <- sapply(all.clusters, function(x){as.numeric(x$z[,1])})
z2 <- sapply(all.clusters, function(x){as.numeric(x$z[,2])})
calls <- sapply(all.clusters, function(x){as.numeric(x$classification)})

# Look at the data in log likelihood space

hist(log10(z1/z2), breaks=10000, xlim=c(-5, 5))

# Bad calls are likely in the middle around 0 from -1 to 1. Turn them into NA.

indNA <- which(log10(z1/z2)>-1&log10(z1/z2)<1, arr.ind=T)

calls[indNA] <- NA

# Create SNPs x Strains matrix

rownames(calls) <- as.character(unique(gg$strain))

calls <- t(calls)

# Polarize the genome by N2 (Bristol, "B") or CB4856 (Hawaiian, "H") alleles based on QG1

qg <- as.numeric(calls[,"QG1"])

qg.polarize <- function(col){
    BH <- ifelse(col == qg, "B", "H")
    return(BH)
}

n.calls <- apply(calls, 2, qg.polarize)

# Check that none of the QG1 SNPs were misclassified as "H". The following line should = 0.

sum(n.calls[,"QG1"] == "H", na.rm=T)

# Get the chromosome and position values for each SNP

chrom <- as.numeric(gg$chr[match(rownames(n.calls), gg$snp)])
pos <- gg$pos[match(rownames(n.calls), gg$snp)]

# Order the calls based on chomosome and position

ord.calls <- n.calls[order(chrom, pos),]
ord.chrom <- chrom[order(chrom, pos)]
ord.pos <- pos[order(chrom, pos)]

ord.id <- paste(ord.chrom, ord.pos, rownames(ord.calls), sep="_")

rownames(ord.calls) <- ord.id

# Add in identifying info for the RIAILs

rils <- paste("QX", seq(240, 598), sep="")

ril.ind <- match(rils, colnames(ord.calls))

# Make the final genotype table for the the cross object and save it to be incorporated into a cross object later

forqtl <- data.frame(id = rownames(ord.calls[,ril.ind]), chrom=ord.chrom, ord.calls[,ril.ind], stringsAsFactors=F)

write.table(forqtl, file="GenotypeProcessing/fullgeno_newrils.csvsr", sep="\t", row.names=F, quote=F)

# Make a dummy phenotype file for the cross object

pheno <- forqtl[c(1,2),-2]

pheno[1,-1] <- rep(5, ncol(pheno)-1) 
pheno[2,-1] <- rep(1, ncol(pheno)-1) 

write.table(pheno, file="GenotypeProcessing/pheno_newrils.csvsr", sep="\t", row.names=F, quote=F)

######################## Begin testing the cross object ########################

# Set the location of the phenotype and genotype files

gfile <- "GenotypeProcessing/fullgeno_newrils.csvsr"
pfile <- "GenotypeProcessing/pheno_newrils.csvsr"

# Define a function to read genotype and phenotype info into a cross object for R/qtl

readCross = function(g,p) {
    return( 
        read.cross(format='csvsr', 
                   genfile = g,
                   phefile = p,
                   na.strings='NA', genotypes=c('B','H'), alleles=c('B','H'), 
                   sep='\t', estimate.map=FALSE, comment.char='')) 
}

# Create the cross object with readCross

cr <- readCross(gfile, pfile)

# Plot all of the genotype skews for "BB"

gt <- geno.table(cr, scanone=T)

ggplot(data=gt) + aes(x=pos, y=BB) + geom_point() + facet_grid(chr~.) +
    geom_hline(yintercept=0.5, color="red") + 
    geom_hline(yintercept=0.45, color="blue") + 
    geom_hline(yintercept=0.4, color="orange") + 
    geom_hline(yintercept=0.35, color="purple") + 
    geom_hline(yintercept=0.3, color="green") + 
    geom_hline(yintercept=0.25, color="yellow")

# Make allele frequency cuts individually for each chromosome for "BB"

cr$geno$"1"$data[,colnames(cr$geno$"1"$data) %in% rownames(gt[which(gt$BB<0.35),])] <- NA
cr$geno$"2"$data[,colnames(cr$geno$"2"$data) %in% rownames(gt[which(gt$BB<0.35),])] <- NA
cr$geno$"3"$data[,colnames(cr$geno$"3"$data) %in% rownames(gt[which(gt$BB<0.3),])] <- NA
cr$geno$"4"$data[,colnames(cr$geno$"4"$data) %in% rownames(gt[which(gt$BB<0.35),])] <- NA
cr$geno$"5"$data[,colnames(cr$geno$"5"$data) %in% rownames(gt[which(gt$BB<0.3),])] <- NA
cr$geno$"6"$data[,colnames(cr$geno$"6"$data) %in% rownames(gt[which(gt$BB<0.3),])] <- NA

# Plot all of the genotype skews for "BH"

ggplot(data=gt) + aes(x=pos, y=BH) + geom_point() + facet_grid(chr~.) +
    geom_hline(yintercept=0.5, color="red") + 
    geom_hline(yintercept=0.45, color="blue") + 
    geom_hline(yintercept=0.4, color="orange") + 
    geom_hline(yintercept=0.35, color="purple") + 
    geom_hline(yintercept=0.3, color="green") + 
    geom_hline(yintercept=0.25, color="yellow")

# Make allele frequency cuts individually for each chromosome for "BH"

cr$geno$"2"$data[,colnames(cr$geno$"2"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"3"$data[,colnames(cr$geno$"3"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"4"$data[,colnames(cr$geno$"4"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"5"$data[,colnames(cr$geno$"5"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"6"$data[,colnames(cr$geno$"6"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA

# Eliminate addition troublesome SNPs that imputation does not fix

cr$geno$"1"$data[,colnames(cr$geno$"1"$data) == "1_13432384_UCE1-1433"] <- NA
cr$geno$"5"$data[,colnames(cr$geno$"5"$data) == "5_19393678_UCE5-3001"] <- NA
cr$geno$"5"$data[,colnames(cr$geno$"5"$data) == "5_14503422_CE5-218"] <- NA
cr$geno$"6"$data[,colnames(cr$geno$"6"$data) == "6_1663691_UCE6-686"] <- NA
cr$geno$"6"$data[,colnames(cr$geno$"6"$data) == "6_12615249_UCE6-1285"] <- NA

# Get all of the SNP names

allSnps <- c(colnames(cr$geno$"1"$data), colnames(cr$geno$"2"$data), colnames(cr$geno$"3"$data), colnames(cr$geno$"4"$data), colnames(cr$geno$"5"$data), colnames(cr$geno$"6"$data))

# Get the set of SNPs not in the Rockman set (Rockman, et al. 2009) and drop these SNPs from the cross object

keepSnps <- c()

for(i in rockmanSNPs){
    i <- paste0(i, "$")
    keepSnps <- append(keepSnps, allSnps[which(grepl(i, allSnps))])
}

removeSnps <- allSnps[!(allSnps %in% keepSnps)]

cr <- drop.markers(cr, removeSnps)

gt3 <- geno.table(cr, scanone=T)

# Impute genotypes on the cross object

icr <- argmax.geno(cr, error.prob=.00001)

# Pull out the imputed genotypes along with position information and save the new genotype table

genoFile <- lapply(icr$geno, function(x){x$argmax})
genoFile <- do.call(cbind, genoFile)
rownames(genoFile) <- paste0("QX", 240:598)
genoFile <- t(genoFile)
id <- sapply(rownames(genoFile), function(x){strsplit(x, "_")[[1]][3]})
chrom <- sapply(rownames(genoFile), function(x){strsplit(x, "_")[[1]][1]})
genoFile <- data.frame(lapply(data.frame(genoFile), function(x){ifelse(x==1, "B", "H")}))
genoFile <- data.frame(cbind(id, chrom, genoFile))
write.table(genoFile, file="GenotypeProcessing/imputedgeno_newrils.csvsr", sep="\t", row.names=F, quote=F)

# Make new cross object from imputed genotype data

pfile <- "GenotypeProcessing/pheno_newrils.csvsr"
gfile <- "GenotypeProcessing/imputedgeno_newrils.csvsr"

cr2 <- readCross(gfile, pfile)

cr.map <- est.map(cr2, verbose=T, error.prob=0.06, tol=1e-5)

plot.map(cr.map)

# Pull out the genotype table

gt2 <- geno.table(cr2, scanone=T)

# Plot the genotype skews. No spikes in the data indicate issues were fixed by imputation.

ggplot(data=gt2) + aes(x=pos, y=BB) + geom_line() + geom_point() +
    facet_grid(chr~., scales="free_x") +
    xlab("Genomic position (Mb)") + ylab("Bristol allele frequency") +
    theme(plot.title=element_text(size=18, face="bold"), 
          axis.title.x = element_text(vjust=0, size=16, face="bold"), 
          legend.position="none", 
          axis.title.y=element_text(size=16, angle=90, face="bold"), 
          axis.text.x=element_text(size=14, face="bold", color="black"), 
          axis.text.y=element_text(size=14, face="bold", color="black"), 
          strip.text.x=element_text(size=14, face="bold"),
          strip.text.y=element_text(size=14, face="bold", angle=90))

# Rename the cross object to work with the mapping pipeline

N2xCB4856.cross <- cr2

save(N2xCB4856.cross, file="Mapping/N2xCB4856_RIAILs_Rqtlfiles.RData")
