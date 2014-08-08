library(plyr)
library(dplyr)
library(ggplot2)
#Read in the table with matt rockmans SNP set
rockman <- read.table("~/Downloads/RILS.txt", header=TRUE)
rockmanSNPs <- rockman$SNP_ID
troubleSNPs <- c("UCE6-686", "UCE4-528", "CE4-241", "UCE6-712", "UCE6-1285", "CE5-218")

#2012 12 07 I will process my RIAILs using EM and clean up with rQTL

load("~/Downloads/rawplates.RData")

#Now, I need to put in strain assignments for the raw data from Rockefeller.

strains <- read.csv("~/Downloads/GG_strains.csv", header=F, stringsAsFactors=F)[,1]

plate1$strains <- rep(strains[1:96], each=1536)
plate2$strains <- rep(strains[97:192], each=1536)
plate3$strains <- rep(strains[289:384], each=1536)
plate4$strains <- rep(strains[193:288], each=1536)

#Note, plate 3 and plate 4 assignments are switched above to fix the plate mixup at Rockefeller. I can only
#really see mixups with plate 4 because there are non-RIAILs on that plate. Plates 1 and 2 are correct.

plate1$plate <- 1
plate2$plate <- 2
plate3$plate <- 3
plate4$plate <- 4

all <- rbind(plate1, plate2, plate3, plate4)

#Reduce to only useful data

gg <- with(all, data.frame(strain = strains, snp=SNP.Name, chr=Chr, pos=Position, x=X.Raw, y=Y.Raw))

#Here is where you can look at the plots of raw values to see how genotypes were assigned. 

gg %>% filter(., snp == "CE3-140") %>% ggplot(.) + aes(x=x, y=y) + geom_point()
gg %>% ggplot(.) + aes(x=x, y=y) + geom_point()

#Genotypes at each SNP can be assigned as a mixture of two normal distributions. Individuals outside of defined groups
#can be called NA and cleaned up with the HMM later.

library(mclust)

#Make an error vector for the tryCatch

mclust.error <- Mclust(data=gg[gg$snp=="CE1-28", c("x","y")], G=2, modelNames=c("VVV"), warn=T)

mclust.error$classification[1:384] <- rep(NA, 384)

#Assign to genotypes using EM algorithm. It takes a few minutes.

all.clusters <- dlply(.data=gg, .variables=.(snp), .fun=function(D) {
  tryCatch({out<-Mclust(data=D[,c("x","y")], G=2, modelNames=c("VVV"), warn=T)
            #out$data<-D[,c("strain", "snp", "chr", "pos")]
            return(out)}, error=function(x){print(D$snp[1]); return(mclust.error)})
})

#Bad marker: UCE2-2369

#Turn into numeric

z1<-sapply(all.clusters, function(x){as.numeric(x$z[,1])})
z2<-sapply(all.clusters, function(x){as.numeric(x$z[,2])})
calls<-sapply(all.clusters, function(x){as.numeric(x$classification)})

#Look at the data in log likelihood space

hist(log10(z1/z2), breaks=10000, xlim=c(-5, 5))

#Bad calls are likely in the middle around 0 from -1 to 1. Turn them into NA

sum(log10(z1/z2)>-1&log10(z1/z2)<1)/length(z1)

indNA <- which(log10(z1/z2)>-1&log10(z1/z2)<1, arr.ind=T)

calls[indNA] <- NA

#Now, look to see how genotype calls look.

rownames(calls) <- as.character(unique(gg$strain))

calls <- t(calls)

qg <- as.numeric(calls[,"QG1"])

qg.polarize <- function(col){
  BH <- ifelse(col == qg, "B", "H")
  return(BH)
}

n.calls <- apply(calls, 2, qg.polarize)

sum(n.calls[,"QG1"] == "H", na.rm=T)

chrom <- as.numeric(gg$chr[match(rownames(n.calls), gg$snp)])
pos <- gg$pos[match(rownames(n.calls), gg$snp)]

ord.calls <- n.calls[order(chrom, pos),]

plot(as.numeric(as.factor(ord.calls[,"QX1436 6X N2"])), type="b")

ord.chrom <- chrom[order(chrom, pos)]
ord.pos <- pos[order(chrom, pos)]

ord.id <- paste(ord.chrom, ord.pos, rownames(ord.calls), sep="_")

rownames(ord.calls) <- ord.id

rils <- paste("QX", seq(240, 598), sep="")

ril.ind <- match(rils, colnames(ord.calls))

forqtl <- data.frame(id = rownames(ord.calls), chrom=ord.chrom, ord.calls[,ril.ind], stringsAsFactors=F)

write.table(forqtl, file="~/Downloads/fullgeno_newrils.csvsr", sep="\t", row.names=F, quote=F)

pheno <- forqtl[c(1,2),-2]

pheno[1,-1] <- rep(5, ncol(pheno)-1) 
pheno[2,-1] <- rep(1, ncol(pheno)-1) 

write.table(pheno, file="~/Downloads/pheno_newrils.csvsr", sep="\t", row.names=F, quote=F)

##############------------------------------------------########################

library(qtl)

gfile <- "~/Downloads/fullgeno_newrils.csvsr"
pfile <- "~/Downloads/pheno_newrils.csvsr"

gfile <- "~/Downloads/imputedgeno.csvsr"

readCross = function(g,p) {   return( 
  read.cross(format='csvsr', 
             genfile = g,
             phefile = p,
             na.strings='NA', genotypes=c('B','H'), alleles=c('B','H'), 
             sep='\t', estimate.map=FALSE, comment.char='')   ) 
}

cr <- readCross(gfile, pfile)

#Let's eliminate genotype skews

gt <- geno.table(cr, scanone=T)
head(gt)

ggplot(data=gt) + aes(x=pos, y=BB) + geom_point() + facet_grid(chr~.) +
  geom_hline(yintercept=0.5, color="red") + 
  geom_hline(yintercept=0.45, color="blue") + 
  geom_hline(yintercept=0.4, color="orange") + 
  geom_hline(yintercept=0.35, color="purple") + 
  geom_hline(yintercept=0.3, color="green") + 
  geom_hline(yintercept=0.25, color="yellow")

#Make allele frequency cuts...BB<0.25 is a good general cutoff

cr$geno$"1"$data[,colnames(cr$geno$"1"$data) %in% rownames(gt[which(gt$BB<0.35),])] <- NA
cr$geno$"2"$data[,colnames(cr$geno$"2"$data) %in% rownames(gt[which(gt$BB<0.35),])] <- NA
cr$geno$"3"$data[,colnames(cr$geno$"3"$data) %in% rownames(gt[which(gt$BB<0.3),])] <- NA
cr$geno$"4"$data[,colnames(cr$geno$"4"$data) %in% rownames(gt[which(gt$BB<0.35),])] <- NA
cr$geno$"5"$data[,colnames(cr$geno$"5"$data) %in% rownames(gt[which(gt$BB<0.3),])] <- NA
cr$geno$"6"$data[,colnames(cr$geno$"6"$data) %in% rownames(gt[which(gt$BB<0.3),])] <- NA

gt2 <- geno.table(cr, scanone=T)

ggplot(data=gt2) + aes(x=pos, y=BB) + geom_point() + facet_grid(chr~.) +
  geom_hline(yintercept=0.5, color="red") + 
  geom_hline(yintercept=0.45, color="blue") + 
  geom_hline(yintercept=0.4, color="orange") + 
  geom_hline(yintercept=0.35, color="purple") + 
  geom_hline(yintercept=0.3, color="green") + 
  geom_hline(yintercept=0.25, color="yellow")

ggplot(data=gt2) + aes(x=pos, y=BH) + geom_point() + facet_grid(chr~.) +
  geom_hline(yintercept=0.5, color="red") + 
  geom_hline(yintercept=0.45, color="blue") + 
  geom_hline(yintercept=0.4, color="orange") + 
  geom_hline(yintercept=0.35, color="purple") + 
  geom_hline(yintercept=0.3, color="green") + 
  geom_hline(yintercept=0.25, color="yellow")

cr$geno$"2"$data[,colnames(cr$geno$"2"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"3"$data[,colnames(cr$geno$"3"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"4"$data[,colnames(cr$geno$"4"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"5"$data[,colnames(cr$geno$"5"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA
cr$geno$"6"$data[,colnames(cr$geno$"6"$data) %in% rownames(gt[which(gt$BH<0.4),])] <- NA

cr$geno$"1"$data[,colnames(cr$geno$"1"$data) == "1_13432384_UCE1-1433"] <- NA
cr$geno$"5"$data[,colnames(cr$geno$"5"$data) == "5_19393678_UCE5-3001"] <- NA
cr$geno$"6"$data[,colnames(cr$geno$"6"$data) == "6_1663691_UCE6-686"] <- NA


gt3 <- geno.table(cr, scanone=T)

ggplot(data=gt3) + aes(x=pos, y=BH) + geom_point() + facet_grid(chr~.) +
  geom_hline(yintercept=0.5, color="red") + 
  geom_hline(yintercept=0.45, color="blue") + 
  geom_hline(yintercept=0.4, color="orange") + 
  geom_hline(yintercept=0.35, color="purple") + 
  geom_hline(yintercept=0.3, color="green") + 
  geom_hline(yintercept=0.25, color="yellow")












cr2 <- cr

allSnps <- c(colnames(cr$geno$"1"$data), colnames(cr$geno$"2"$data), colnames(cr$geno$"3"$data), colnames(cr$geno$"4"$data), colnames(cr$geno$"5"$data), colnames(cr$geno$"6"$data))

keepSnps <- c()

for(i in rockmanSNPs){
    i <- paste0(i, "$")
    keepSnps <- append(keepSnps, allSnps[which(grepl(i, allSnps))])
}

removeSnps <- allSnps[!(allSnps %in% keepSnps)]

cr2 <- drop.markers(cr2, removeSnps)

# argmax.geno() here, play with error prob
icr <- argmax.geno(cr2, error.prob=.00001)




View(filter(gt, chr==1))











a <- data.frame(icr$geno$"1"$argmax)


gt4 <- geno.table(icr, scanone=T)
gt4$chr <- factor(gt4$chr, labels=c("I", "II", "III", "IV", "V", "X"))

gt4$SNP <- do.call(rbind, strsplit(rownames(gt4), "_"))[,3]

gt4 <- filter(gt4, SNP %in% rockmanSNPs)

load("markers.Rda")

markers2 <- select(markers, -chr.num, -chr.roman)
gt4 <- merge(gt4, markers2, by="SNP")

write.table(gt4, file="~/Downloads/newrils_genotable.txt", sep="\t", quote=F, row.names=F, col.names=T)

ggplot(data=gt4) + aes(x=WS185.pos/1e6, y=BB) + geom_line() + geom_point() +
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

#Skews are effectively gone.

#Now, check for similarities in genotype

cg <- comparegeno(cr)

hist(cg, breaks=100, xlab="Fraction of genome shared", main="")
abline(v=0.5, col="red", lty=2, lwd=3)
rug(cg)

which(cg>0.9, arr.ind=T)

#There are no similar individuals.

#Now, check recombination fraction.

cr <- est.rf(cr)

checkAlleles(cr) #1_316117_UCE1-530 is likely switched.

g <- pull.geno(cr, 1)
g[,"1_316117_UCE1-530"] <- abs(g[,"1_316117_UCE1-530"] - 3)
cr$geno[[1]]$data <- g

#This has to be done again for no apparent errors to be found
cr <- est.rf(cr)

checkAlleles(cr) #No apparent problems. It's fixed.

plotRF(cr)

#There are some markers out of place, but better than Matt's set.

cr.map <- est.map(cr, verbose=T, error.prob=0.06, tol=1e-5)

plot.map(cr.map)

#Chromosome map looks pretty bad. Let's try to impute and see what happens.

imputeMissingGeno  = function(cross) {
  impcross = argmax.geno(cross, error.prob=0.06)
  impcross$geno = lapply(impcross$geno, function(x){x$data = x$argmax; x$argmax=NULL; return(x) })
  #dm = findDupMarkers(impcross, exact.only = FALSE , adjacent.only = FALSE )
  #names(dm)=NULL
  #impcross = drop.markers(impcross, unlist(dm))
  return(impcross) }

icr <- imputeMissingGeno(cr)

icr <- est.rf(icr)

checkAlleles(icr) #No apparent problems.

plotRF(icr) #There might be some more LD between 1 and 2, but otherwise way better.

icr.map <- est.map(icr, verbose=T, error.prob=0.06, tol=1e-5)

newrils.impute <- icr

plot.map(icr.map) #Map looks great. Imputation smoothed over bad SNPs from before.

#I need to make a position file for distribution to labs to do their own mapping.
library(stringr)

snplist <- c(colnames(icr$geno$"1"$data), colnames(icr$geno$"2"$data), colnames(icr$geno$"3"$data),
      colnames(icr$geno$"4"$data), colnames(icr$geno$"5"$data), colnames(icr$geno$"6"$data))

sp.snplist <- str_split_fixed(snplist, "_", 3)

sp.snplist[,1] <- ifelse(sp.snplist[,1]=="6", "X", sp.snplist[,1])

snps <- data.frame(fullsnp = as.character(snplist), snp = as.character(sp.snplist[,3]), chr = as.factor(sp.snplist[,1]), pos = as.numeric(sp.snplist[,2]), stringsAsFactors=F)

setwd("~/Dropbox/GoldenGate for new RIAILs/")

save(snps, newrils.impute, file="newrils.RData")

#Let's try a mapping

names(newrils.impute$geno) <- chrnames

newrils.impute <- calc.genoprob(newrils.impute)

newqx <- paste("QX", seq(240, 598), sep="")

red.newrils <- subset(newrils.impute, chr=chrnames, id = newqx)

red.pq <- subset(all.pq, chr=chrnames, ind = newqx)

comb.new <- red.newrils

comb.new$pheno <- red.pq$pheno

trait <- scanone(comb.new, pheno.col=5, model="np")
plot(trait)






genoFile <- lapply(icr$geno, function(x){x$argmax})
genoFile <- do.call(cbind, genoFile)
rownames(genoFile) <- paste0("QX", 240:598)
genoFile <- t(genoFile)
id <- sapply(rownames(genoFile), function(x){strsplit(x, "_")[[1]][3]})
chrom <- sapply(rownames(genoFile), function(x){strsplit(x, "_")[[1]][1]})
genoFile <- data.frame(lapply(data.frame(genoFile), function(x){ifelse(x==1, "B", "H")}))
genoFile <- data.frame(cbind(id, chrom, genoFile))
write.table(genoFile, file="~/Downloads/imputedgeno.csvsr", sep="\t", row.names=F, quote=F)
