########################################################################################
########################################################################################
########################################################################################
########################################################################################
# For RIAILs0 figures, here is the plan. All figures should be in theme_bw() to be compatible with publication.
# 
# Figure 1: Allele frequencies across each of the six chromosomes with a dotted red line at a frequency of 0.5. Matt's frequencies can be shown in light gray. Our frequencies can be shown in black. -----DONE-----
# 
# Figure 2: Power calculations for larger 359 RIAIL collection. -----DONE-----
# 
# Figure 3: Pictorial overview of the HTA pipeline -----DONE-----
# 
# Figure 4a: Boxplots for norm.n of N2, CB4856, and RIAILs in control conditions. N2 in orange, CB in blue, and RIAILs in gray. -----DONE----- 
# Figure 4b: Mapping of norm.n in control conditions. Line should be weight 1. X-axis is physical position, free_x on faceting by chromosome, y-axis is LOD score, peaks are shown as upside down red triangles -----DONE-----
# 
# I'm imagining a and b will be on top, c and d in middle, and e and f on bottom for a nearly full page figure.
# Figure 5a: Histogram of median TOF in control conditions -----DONE-----
# Figure 5b: Map of median TOF -----DONE-----
# Figure 5c: Histogram of median EXT in control conditions -----DONE-----
# Figure 5d: Map of median EXT -----DONE-----
# Figure 5e: Histogram of median norm.EXT in control conditions -----DONE-----
# Figure 5f: Map of median norm.EXT -----DONE-----
# 
# Here is our cool figure as discussed before. Each quantile or IQR or trait will be its own color on the central 96-well trait array. 
# Figure 6a: 96-well layout of median norm.EXT histograms with color lines and ranges
# Figure 6b-g: q10, q25, median, q75, q90, IQR or var
# 
# Figure 7: Paraquat dose response
# 
# Figure 8a: histogram of median norm.EXT in control (black) and paraquat (red) conditions
# Figure 8b: Map of paraquat trait that maps
# 
# Figure 9: NIL results (if we still have overlap with real QTL)
# 
# Table 1: QX1430 vs. N2 for genetic incompatibility
# 
# Table 2: Mapping results: Treatment, trait, QTL#, marker, peak phys. position, L CI, R CI, variance explained 
########################################################################################
########################################################################################
########################################################################################
########################################################################################

library(dplyr) #Version 0.2
library(ggplot2)

load("PhenotypeProcessing/PresentationStyle.RData")

## Figure 4A
fig4DataA <- read.csv("Data/MappingPhenotypes.csv")
nData <- fig4DataA %>% filter(drug=="control") %>% select(strain, n) %>% mutate(group=ifelse(strain=="N2", "N2", ifelse(strain=="CB4856", "CB4856", "RIAILs")))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_broodsize/FecundityBoxplot.png", height=3, width=5)

ggplot(nData, aes(x=factor(group), y=n)) + geom_boxplot(aes(fill=factor(group))) + theme_bw() + xlab(NULL) + ylab("Number of offspring") + theme(legend.position="none") + ggtitle("Lifetime Fecundity") + scale_x_discrete(limits=c("N2", "CB4856", "RIAILs"), labels=c("Bristol", "Hawaii", "RIAILs")) + scale_fill_manual(values = c("N2" = "blue","CB4856" = "orange","RIAILs" = "gray")) + presentation

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_broodsize/FecundityMap.png", height=3, width=5)

## Figure 4B
fig4DataBControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.n")
fig4DataBControlMap$chr <- ifelse(fig4DataBControlMap$chr==1, "I",
                                  ifelse(fig4DataBControlMap$chr==2, "II",
                                         ifelse(fig4DataBControlMap$chr==3, "III",
                                                ifelse(fig4DataBControlMap$chr==4, "IV",
                                                       ifelse(fig4DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig4DataBControlMap[!is.na(fig4DataBControlMap$var.exp) & as.character(fig4DataBControlMap$trait)=="control.n",]
ggplot(fig4DataBControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("Lifetime Fecundity in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_broodsize/FecundityMap.png", height=5, width=10)

## Figure 5
fig5Data <- read.csv("Data/MappingPhenotypes.csv")

# A
ggplot(fig5Data[fig5Data$drug=="control",], aes(x = median.TOF)) + geom_bar() + theme_bw() + presentation + ggtitle("Distribution of Median Time of Flight Values\nin Control Conditions") + xlab("Median Time of Flight") + ylab("Count")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/TOFHist.png", height=5, width=10)

# B
fig5DataBControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.median.TOF")
fig5DataBControlMap$chr <- ifelse(fig5DataBControlMap$chr==1, "I",
                                  ifelse(fig5DataBControlMap$chr==2, "II",
                                         ifelse(fig5DataBControlMap$chr==3, "III",
                                                ifelse(fig5DataBControlMap$chr==4, "IV",
                                                       ifelse(fig5DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig5DataBControlMap[!is.na(fig5DataBControlMap$var.exp) & as.character(fig5DataBControlMap$trait)=="control.median.TOF",]
ggplot(fig5DataBControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("Median Time of Flight in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/TOFMap.png", height=5, width=10)


# C
ggplot(fig5Data[fig5Data$drug=="control",], aes(x = median.EXT)) + geom_bar() + theme_bw() + presentation + ggtitle("Distribution of Median Extinction Values\nin Control Conditions") + xlab("Median Extinction") + ylab("Count")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/EXTHist.png", height=5, width=10)

# D
fig5DataDControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.median.EXT")
fig5DataDControlMap$chr <- ifelse(fig5DataDControlMap$chr==1, "I",
                                  ifelse(fig5DataDControlMap$chr==2, "II",
                                         ifelse(fig5DataDControlMap$chr==3, "III",
                                                ifelse(fig5DataDControlMap$chr==4, "IV",
                                                       ifelse(fig5DataDControlMap$chr==5, "V", "X")))))
peaksDF <- fig5DataDControlMap[!is.na(fig5DataDControlMap$var.exp) & as.character(fig5DataDControlMap$trait)=="control.median.EXT",]
ggplot(fig5DataDControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("Median Extinction in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/EXTMap.png", height=5, width=10)


# E
ggplot(fig5Data[fig5Data$drug=="control",], aes(x = median.norm.EXT)) + geom_bar() + theme_bw() + presentation + ggtitle("Distribution of Median Normalized Extinction Values\nin Control Conditions") + xlab("Median Normalized Extinction") + ylab("Count")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/normEXTHist.png", height=5, width=10)

# F
fig5DataFControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.median.norm.EXT")
fig5DataFControlMap$chr <- ifelse(fig5DataFControlMap$chr==1, "I",
                                  ifelse(fig5DataFControlMap$chr==2, "II",
                                         ifelse(fig5DataFControlMap$chr==3, "III",
                                                ifelse(fig5DataFControlMap$chr==4, "IV",
                                                       ifelse(fig5DataFControlMap$chr==5, "V", "X")))))
peaksDF <- fig5DataFControlMap[!is.na(fig5DataFControlMap$var.exp) & as.character(fig5DataFControlMap$trait)=="control.median.norm.EXT",]
ggplot(fig5DataFControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("Median Normalized Extinction in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/normEXTMap.png", height=5, width=10)




## Figure 6
load("Data/RawScoreData.Rda")
summarizedScoreData <- read.csv("Data/ProcessedPhenotypes.csv") %>% select(assay, plate, row, col, strain)
histsData <- rawScoreData %>%
    filter(drug=="control") %>%
    select(assay, plate, row, col, TOF) %>%
    filter(as.numeric(as.character(col)) %in% c(1, 5, 9))
histsData$strain <- sapply(1:nrow(histsData), function(i){
    print(i)
    as.character(summarizedScoreData[summarizedScoreData$assay==histsData$assay[i] & as.numeric(as.character(summarizedScoreData$plate))==as.numeric(as.character(histsData$plate[i])) & summarizedScoreData$row==histsData$row[i] & summarizedScoreData$col==histsData$col[i] , "strain"])
})

histsData2 <- histsData[sapply(1:nrow(histsData), function(x)as.numeric(strsplit(histsData$strain[x], "X")[[1]][2]) > 239),]

set.seed(0)

strains <- data.frame(table(histsData2$strain)) %>% 
    arrange(desc(Freq)) %>% .[sample(1:338, 12), 1] %>%
    as.character(.)
histsData3 <- histsData2[histsData2$strain %in% strains,]
quantileData <- histsData3 %>%
    group_by(strain) %>%
    summarize(q10=quantile(TOF, probs=.10),
           q25=quantile(TOF, probs=.25),
           q50=quantile(TOF, probs=.50),
           mean=mean(TOF),
           q75=quantile(TOF, probs=.75),
           q90=quantile(TOF, probs=.90))
    
ggplot(histsData3, aes(x = TOF)) +
    theme_bw() +
    presentation +
    geom_bar() +
    facet_wrap(~strain, ncol=6) +
    geom_vline(data=quantileData, aes(xintercept=q10), colour = "red") +
    geom_vline(data=quantileData, aes(xintercept=q25), colour = "orange") +
    geom_vline(data=quantileData, aes(xintercept=q50), colour = "yellow") +
#     geom_vline(data=quantileData, aes(xintercept=mean), colour = "green") +
    geom_vline(data=quantileData, aes(xintercept=q75), colour = "blue") +
    geom_vline(data=quantileData, aes(xintercept=q90), colour = "purple") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    xlab("Time of Flight Value") +
    ylab("Count")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/distributions.png", height=7, width=10)

# B
fig6DataBControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.q10.TOF")
fig6DataBControlMap$chr <- ifelse(fig6DataBControlMap$chr==1, "I",
                                  ifelse(fig6DataBControlMap$chr==2, "II",
                                         ifelse(fig6DataBControlMap$chr==3, "III",
                                                ifelse(fig6DataBControlMap$chr==4, "IV",
                                                       ifelse(fig6DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataBControlMap[!is.na(fig6DataBControlMap$var.exp) & as.character(fig6DataBControlMap$trait)=="control.q10.TOF",]
ggplot(fig6DataBControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("10th Quantile of Time of Flight in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q10Map.png", height=5, width=10)

# C
fig6DataCControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.q25.TOF")
fig6DataCControlMap$chr <- ifelse(fig6DataCControlMap$chr==1, "I",
                                  ifelse(fig6DataCControlMap$chr==2, "II",
                                         ifelse(fig6DataCControlMap$chr==3, "III",
                                                ifelse(fig6DataCControlMap$chr==4, "IV",
                                                       ifelse(fig6DataCControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataCControlMap[!is.na(fig6DataCControlMap$var.exp) & as.character(fig6DataCControlMap$trait)=="control.q25.TOF",]
ggplot(fig6DataCControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("25th Quantile of Time of Flight in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q25Map.png", height=5, width=10)

# D
fig6DataDControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.median.TOF")
fig6DataDControlMap$chr <- ifelse(fig6DataDControlMap$chr==1, "I",
                                  ifelse(fig6DataDControlMap$chr==2, "II",
                                         ifelse(fig6DataDControlMap$chr==3, "III",
                                                ifelse(fig6DataDControlMap$chr==4, "IV",
                                                       ifelse(fig6DataDControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataDControlMap[!is.na(fig6DataDControlMap$var.exp) & as.character(fig6DataDControlMap$trait)=="control.median.TOF",]
ggplot(fig6DataDControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("Median Time of Flight in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/medianMap.png", height=5, width=10)

# E
fig6DataEControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.q75.TOF")
fig6DataEControlMap$chr <- ifelse(fig6DataEControlMap$chr==1, "I",
                                  ifelse(fig6DataEControlMap$chr==2, "II",
                                         ifelse(fig6DataEControlMap$chr==3, "III",
                                                ifelse(fig6DataEControlMap$chr==4, "IV",
                                                       ifelse(fig6DataEControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataEControlMap[!is.na(fig6DataEControlMap$var.exp) & as.character(fig6DataEControlMap$trait)=="control.q75.TOF",]
ggplot(fig6DataEControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("75th Quantile of Time of Flight in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q75Map.png", height=5, width=10)

# F
fig6DataFControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.q90.TOF")
fig6DataFControlMap$chr <- ifelse(fig6DataFControlMap$chr==1, "I",
                                  ifelse(fig6DataFControlMap$chr==2, "II",
                                         ifelse(fig6DataFControlMap$chr==3, "III",
                                                ifelse(fig6DataFControlMap$chr==4, "IV",
                                                       ifelse(fig6DataFControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataFControlMap[!is.na(fig6DataFControlMap$var.exp) & as.character(fig6DataFControlMap$trait)=="control.q90.TOF",]
ggplot(fig6DataFControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("90th Quantile of Time of Flight in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q90Map.png", height=5, width=10)

# G
fig6DataGControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="control.iqr.TOF")
fig6DataGControlMap$chr <- ifelse(fig6DataGControlMap$chr==1, "I",
                                  ifelse(fig6DataGControlMap$chr==2, "II",
                                         ifelse(fig6DataGControlMap$chr==3, "III",
                                                ifelse(fig6DataGControlMap$chr==4, "IV",
                                                       ifelse(fig6DataGControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataBControlMap[!is.na(fig6DataBControlMap$var.exp) & as.character(fig6DataBControlMap$trait)=="control.iqr.TOF",]
ggplot(fig6DataBControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("Interquartile Range of Time of Flight in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/varMap.png", height=5, width=10)



## Figure 8
# A
load("Data/RawScoreData.Rda")
ggplot(rawScoreData, aes(x = norm.EXT)) + geom_bar(aes(fill = drug), position="dodge") + scale_fill_manual(values=c("control"="black", "paraquat"="red")) + presentation + xlab("Normalized Extinction Value") + ylab("Count")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normEXT/normEXThist.png", height=5, width=7)

# B
fig8DataBControlMap <- read.csv("Mapping/MappingResults.csv") %>% filter(trait=="paraquat.mean.norm.EXT")
fig8DataBControlMap$chr <- ifelse(fig8DataBControlMap$chr==1, "I",
                                  ifelse(fig8DataBControlMap$chr==2, "II",
                                         ifelse(fig8DataBControlMap$chr==3, "III",
                                                ifelse(fig8DataBControlMap$chr==4, "IV",
                                                       ifelse(fig8DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig8DataBControlMap[!is.na(fig8DataBControlMap$var.exp) & as.character(fig8DataBControlMap$trait)=="paraquat.mean.norm.EXT",]
ggplot(fig8DataBControlMap) +
    theme_bw() +
    presentation +
    geom_line(aes(x=pos/1e6, y=LOD), size=1) +
    facet_grid(.~chr, scales="free_x") +
    xlab("Position (Mb)") +
    ylab("LOD") +
    ggtitle("Mean Normalized Extinction in Control Conditions") +
    geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
    geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
    geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
    geom_text(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.25*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normEXT/PQnormEXTmap.png", height=5, width=10)