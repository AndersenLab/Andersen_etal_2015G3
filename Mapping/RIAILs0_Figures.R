########################################################################################
########################################################################################
########################################################################################
########################################################################################
# Figures for Andersen et al. 2015 RIAILs0 manuscript 
########################################################################################
########################################################################################
########################################################################################
########################################################################################

library(dplyr) #Version 0.4
library(ggplot2)

## Figure 1, overview of HTA assay, made separately

## Figure 2, Fecundity in control conditions and mapping
# A
fig4DataA <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/ProcessedPhenotypes.csv")
nData <- fig4DataA %>% filter(drug=="control") %>% select(strain, resid.a.n) %>% mutate(group=ifelse(strain=="N2", "N2", ifelse(strain=="CB4856", "CB4856", "RIAILs")))

ggplot(nData, aes(x=factor(group), y=resid.a.n)) + 
  geom_boxplot(outlier.size=1, aes(fill=factor(group))) + 
  theme_bw() + 
  labs(x="", y="Normalized Fecundity", title="") +
  theme(legend.position="none") + 
  scale_x_discrete(limits=c("N2", "CB4856", "RIAILs"), labels=c("Bristol", "CB4856", "RIAILs")) + 
  scale_fill_manual(values = c("N2" = "orange","CB4856" = "blue","RIAILs" = "gray")) + 
  theme(axis.text.x = element_text(size=12, face="bold", color="black"),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_broodsize/FecundityBoxplot.tiff", height=3, width=3, units="in", dpi=300)

## B
fig4DataBControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.n")
#fig4DataBControlMap$chr <- ifelse(fig4DataBControlMap$chr==1, "I",
                                 # ifelse(fig4DataBControlMap$chr==2, "II",
                                  #       ifelse(fig4DataBControlMap$chr==3, "III",
                                   #             ifelse(fig4DataBControlMap$chr==4, "IV",
                                    #                   ifelse(fig4DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig4DataBControlMap[!is.na(fig4DataBControlMap$var.exp) & as.character(fig4DataBControlMap$trait)=="control.resid.a.n",]

ggplot(fig4DataBControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
  geom_text(data=peaksDF, size=3, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%"))) +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_broodsize/FecundityMap.tiff", height=3, width=4.5, units="in", dpi=300)

###----------------------------------------------------------------------------###
## Figure 3, because this code randomly selects 10 RIAILs to make the histograms. Do not re-run unless you want to redo the Figure.

load("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/RawScoreData.Rda")

summarizedScoreData <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/ProcessedPhenotypes.csv") %>% select(assay, plate, row, col, strain)

histsData <- rawScoreData %>%
  filter(drug=="control") %>%
  select(assay, plate, row, col, TOF) %>%
  filter(as.numeric(as.character(col)) %in% c(1, 5, 9))

histsData$assay <- as.character(histsData$assay)
histsData$plate <- as.integer(histsData$plate)
histsData$row <- as.character(histsData$row)
histsData$col <- as.integer(histsData$col)

histsData <- histsData %>% mutate(ind = paste(assay, plate, row, col, sep=""))
summarizedScoreData <- summarizedScoreData %>% mutate(ind = paste(assay, plate, row, col, sep=""))

histData2 <- merge(histsData, summarizedScoreData, by="ind") %>% select(assay.x, plate.x, row.x, col.x, TOF, strain) %>% 
  rename(assay = assay.x, plate = plate.x, row = row.x, col = col.x) %>% mutate(ind = as.numeric(str_split_fixed(strain, "QX", 2)[,2])) %>%
  filter(ind > 239) %>% select(-ind)

set.seed(0)

strains <- data.frame(table(histData2$strain)) %>% 
  arrange(desc(Freq)) %>% .[sample(1:338, 10), 1] %>%
  as.character(.)
histsData3 <- histData2[histData2$strain %in% strains,]
quantileData <- histsData3 %>%
  group_by(strain) %>%
  summarize(q10=quantile(TOF, probs=.10),
            q25=quantile(TOF, probs=.25),
            q50=quantile(TOF, probs=.50),
            mean=mean(TOF),
            q75=quantile(TOF, probs=.75),
            q90=quantile(TOF, probs=.90))

ggplot(histsData3, aes(x = TOF)) +
  geom_bar() +
  geom_vline(data=quantileData, aes(xintercept=q10), colour = "red") +
  geom_vline(data=quantileData, aes(xintercept=q25), colour = "orange") +
  geom_vline(data=quantileData, aes(xintercept=q50), colour = "yellow") +
  #     geom_vline(data=quantileData, aes(xintercept=mean), colour = "green") +
  geom_vline(data=quantileData, aes(xintercept=q75), colour = "blue") +
  geom_vline(data=quantileData, aes(xintercept=q90), colour = "purple") +
  facet_wrap(~strain, ncol=5) +
  xlim(0,750) +
  labs(x="Body length (µm)", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, color="black", angle=45, hjust=1),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/distributions.tiff", height=3, width=9, units="in", dpi=300)

# B
fig6DataBControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.q10.TOF")
# fig6DataBControlMap$chr <- ifelse(fig6DataBControlMap$chr==1, "I",
#                                   ifelse(fig6DataBControlMap$chr==2, "II",
#                                          ifelse(fig6DataBControlMap$chr==3, "III",
#                                                 ifelse(fig6DataBControlMap$chr==4, "IV",
#                                                        ifelse(fig6DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataBControlMap[!is.na(fig6DataBControlMap$var.exp) & as.character(fig6DataBControlMap$trait)=="control.resid.a.q10.TOF",]

ggplot(fig6DataBControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=2) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 1.5, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10, face="bold", color="black"),
        axis.title.y = element_text(size=10, face="bold", color="black"),
        strip.text.x = element_text(size=10, face="bold", color="white"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        strip.background = element_rect(fill = "red"),
        plot.title = element_text(size=10, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q10Map.tiff", height=2, width=3, units="in", dpi=300)

# C
fig6DataCControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.q25.TOF")
# fig6DataCControlMap$chr <- ifelse(fig6DataCControlMap$chr==1, "I",
#                                   ifelse(fig6DataCControlMap$chr==2, "II",
#                                          ifelse(fig6DataCControlMap$chr==3, "III",
#                                                 ifelse(fig6DataCControlMap$chr==4, "IV",
#                                                        ifelse(fig6DataCControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataCControlMap[!is.na(fig6DataCControlMap$var.exp) & as.character(fig6DataCControlMap$trait)=="control.resid.a.q25.TOF",]

ggplot(fig6DataCControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=2) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 1.5, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10, face="bold", color="black"),
        axis.title.y = element_text(size=10, face="bold", color="black"),
        strip.text.x = element_text(size=10, face="bold", color="white"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        strip.background = element_rect(fill = "orange"),
        plot.title = element_text(size=10, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q25Map.tiff", height=2, width=3, units="in", dpi=300)

# D
fig6DataDControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.median.TOF")
# fig6DataDControlMap$chr <- ifelse(fig6DataDControlMap$chr==1, "I",
#                                   ifelse(fig6DataDControlMap$chr==2, "II",
#                                          ifelse(fig6DataDControlMap$chr==3, "III",
#                                                 ifelse(fig6DataDControlMap$chr==4, "IV",
#                                                        ifelse(fig6DataDControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataDControlMap[!is.na(fig6DataDControlMap$var.exp) & as.character(fig6DataDControlMap$trait)=="control.resid.a.median.TOF",]

ggplot(fig6DataDControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=2) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 1.5, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10, face="bold", color="black"),
        axis.title.y = element_text(size=10, face="bold", color="black"),
        strip.text.x = element_text(size=10, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        strip.background = element_rect(fill = "yellow"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/medianMap.tiff", height=2, width=3, units="in", dpi=300)

# E
fig6DataEControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.q75.TOF")
# fig6DataEControlMap$chr <- ifelse(fig6DataEControlMap$chr==1, "I",
#                                   ifelse(fig6DataEControlMap$chr==2, "II",
#                                          ifelse(fig6DataEControlMap$chr==3, "III",
#                                                 ifelse(fig6DataEControlMap$chr==4, "IV",
#                                                        ifelse(fig6DataEControlMap$chr==5, "V", "X")))))
peaksDF <- fig6DataEControlMap[!is.na(fig6DataEControlMap$var.exp) & as.character(fig6DataEControlMap$trait)=="control.resid.a.q75.TOF",]

ggplot(fig6DataEControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=2) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 1.5, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10, face="bold", color="black"),
        axis.title.y = element_text(size=10, face="bold", color="black"),
        strip.text.x = element_text(size=10, face="bold", color="white"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        strip.background = element_rect(fill = "blue"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q75Map.tiff", height=2, width=3, units="in", dpi=300)

# F
fig6DataFControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.q90.TOF")

peaksDF <- fig6DataFControlMap[!is.na(fig6DataFControlMap$var.exp) & as.character(fig6DataFControlMap$trait)=="control.resid.a.q90.TOF",]

ggplot(fig6DataFControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=2) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 1.5, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10, face="bold", color="black"),
        axis.title.y = element_text(size=10, face="bold", color="black"),
        strip.text.x = element_text(size=10, face="bold", color="white"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        strip.background = element_rect(fill = "purple"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/q90Map.tiff", height=2, width=3, units="in", dpi=300)

# G
fig6DataGControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.var.TOF")

peaksDF <- fig6DataGControlMap[!is.na(fig6DataGControlMap$var.exp) & as.character(fig6DataGControlMap$trait)=="control.resid.a.var.TOF",]

ggplot(fig6DataGControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=2) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 1.5, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10, face="bold", color="black"),
        axis.title.y = element_text(size=10, face="bold", color="black"),
        strip.text.x = element_text(size=10, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_TOFdistribution/varMap.tiff", height=2, width=3, units="in", dpi=300)

#####------------------------------------------------------------------------#####
## Figure 4
df <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/ProcessedPhenotypes.csv")

# A
ggplot(df[df$drug=="control",], aes(x = resid.a.median.TOF)) + 
  geom_histogram() + 
  labs(x="Median body length (µm)", y="Count") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/TOFHist.tiff", height=3, width=3, units="in", dpi=300)

# B
fig5DataBControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.median.TOF")
# fig5DataBControlMap$chr <- ifelse(fig5DataBControlMap$chr==1, "I",
#                                   ifelse(fig5DataBControlMap$chr==2, "II",
#                                          ifelse(fig5DataBControlMap$chr==3, "III",
#                                                 ifelse(fig5DataBControlMap$chr==4, "IV",
#                                                        ifelse(fig5DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig5DataBControlMap[!is.na(fig5DataBControlMap$var.exp) & as.character(fig5DataBControlMap$trait)=="control.resid.a.median.TOF",]

ggplot(fig5DataBControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 3, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/TOFMap.tiff", height=3, width=4.5, units="in", dpi=300)

# C
ggplot(df[df$drug=="control",], aes(x = resid.a.median.EXT)) + 
  geom_histogram() + 
  labs(x="Median optical density", y="Count") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/EXTHist.tiff", height=3, width=3, units="in", dpi=300)

# D
fig5DataDControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.median.EXT")
# fig5DataDControlMap$chr <- ifelse(fig5DataDControlMap$chr==1, "I",
#                                   ifelse(fig5DataDControlMap$chr==2, "II",
#                                          ifelse(fig5DataDControlMap$chr==3, "III",
#                                                 ifelse(fig5DataDControlMap$chr==4, "IV",
#                                                        ifelse(fig5DataDControlMap$chr==5, "V", "X")))))
peaksDF <- fig5DataDControlMap[!is.na(fig5DataDControlMap$var.exp) & as.character(fig5DataDControlMap$trait)=="control.resid.a.median.EXT",]

ggplot(fig5DataDControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 3, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  scale_y_continuous(breaks = c(1:4*2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/EXTMap.tiff", height=3, width=4.5, units="in", dpi=300)

# E

ggplot(df[df$drug=="control",], aes(x = resid.a.median.norm.EXT)) + 
  geom_histogram() + 
  labs(x="Median normalized opt. density", y="Count") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/normEXTHist.tiff", height=3, width=3, units="in", dpi=300)

# F
fig5DataFControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="control.resid.a.median.norm.EXT")
# fig5DataFControlMap$chr <- ifelse(fig5DataFControlMap$chr==1, "I",
#                                   ifelse(fig5DataFControlMap$chr==2, "II",
#                                          ifelse(fig5DataFControlMap$chr==3, "III",
#                                                 ifelse(fig5DataFControlMap$chr==4, "IV",
#                                                        ifelse(fig5DataFControlMap$chr==5, "V", "X")))))
peaksDF <- fig5DataFControlMap[!is.na(fig5DataFControlMap$var.exp) & as.character(fig5DataFControlMap$trait)=="control.resid.a.median.norm.EXT",]

ggplot(fig5DataFControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_text(data=peaksDF, size= 3, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%")), fontface="bold") +
  facet_grid(.~chr, scales="free_x") +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normalizedEXT/normEXTMap.tiff", height=3, width=4.5, units="in", dpi=300)

#####------------------------------------------------------------------------#####
## Figure 5

load(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_doseresponses/oldpqdose.RData")

ggplot(proc.pq3) + aes(x=conc, y=mean.n, color=strain) + geom_line() +
  labs(x="Paraquat concentration (mM)", y="Mean fecundity") +
  scale_color_manual(values=c("orange", "blue", "red", "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=0))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_doseresponses/pqdose_n.tiff", height=3, width=6, units="in", dpi=300)

#####------------------------------------------------------------------------#####
## Figure 6

# A
load("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/RawScoreData.Rda")
summarizedScoreData <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/ProcessedPhenotypesFull.csv") %>% select(assay, plate, row, col, strain)

pqcomp <- rawScoreData %>%
  select(assay, plate, drug, row, col, norm.EXT)

pqcomp$assay <- as.character(pqcomp$assay)
pqcomp$plate <- as.integer(pqcomp$plate)
pqcomp$row <- as.character(pqcomp$row)
pqcomp$col <- as.integer(pqcomp$col)

pqcomp2 <- pqcomp %>% mutate(ind = paste(assay, plate, row, col, sep=""))
summarizedScoreData <- summarizedScoreData %>% mutate(ind = paste(assay, plate, row, col, sep=""))

pqcomp3 <- merge(pqcomp2, summarizedScoreData, by="ind") %>% select(assay.x, plate.x, drug, row.x, col.x, norm.EXT, strain) %>% 
  rename(assay = assay.x, plate = plate.x, row = row.x, col = col.x) %>% mutate(ind = as.numeric(str_split_fixed(strain, "QX", 2)[,2])) %>%
  filter(ind > 239) %>% select(-ind)

ggplot(pqcomp3, aes(x = norm.EXT)) + 
  geom_bar(aes(fill = drug), position="dodge") + 
  xlab("Normalized optical density") + ylab("Count") + 
  scale_fill_manual(values=c("control"="black", "paraquat"="red")) + 
  theme_bw() +
  theme(axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=10, face="bold", color="black"),
        axis.title.y = element_text(size=10, face="bold", color="black"),
        strip.text.x = element_text(size=10, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        legend.position="none",
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normEXTpq/normEXThist.tiff", height=3, width=3, units="in", dpi=300)

# B
fig8DataBControlMap <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(trait=="paraquat.resid.mean.norm.EXT")
# fig8DataBControlMap$chr <- ifelse(fig8DataBControlMap$chr==1, "I",
#                                   ifelse(fig8DataBControlMap$chr==2, "II",
#                                          ifelse(fig8DataBControlMap$chr==3, "III",
#                                                 ifelse(fig8DataBControlMap$chr==4, "IV",
#                                                        ifelse(fig8DataBControlMap$chr==5, "V", "X")))))
peaksDF <- fig8DataBControlMap[!is.na(fig8DataBControlMap$var.exp) & as.character(fig8DataBControlMap$trait)=="paraquat.resid.mean.norm.EXT",]

ggplot(fig8DataBControlMap) +
  geom_line(aes(x=pos/1e6, y=LOD), size=0.5) +
  geom_segment(data=peaksDF, aes(x=CI.L.pos/1e6, xend=CI.R.pos/1e6, y=0, yend=0), colour="blue", size=2) +
  geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
  geom_point(data=peaksDF, aes(x=pos/1e6, y = LOD+((1.1*max(LOD)-max(LOD)))), fill="red", shape=25, size=4) +
  geom_text(data=peaksDF, size=3, aes(x=pos/1e6, y = LOD+((1.27*max(LOD)-max(LOD))-(1.1*max(LOD)-max(LOD))), label=paste0(round(100*var.exp, 2), "%"))) +
  labs(x="Position (Mb)", y="LOD") +
  facet_grid(.~chr, scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_normEXTpq/PQresidnormEXTmap.tiff", height=3, width=4.5, units="in", dpi=300)

####----------------------------------------------------#####
## Figure S1
#A - Worms and bubbles

load(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_bubblesROC/bubbleworm.RData")

ggplot(data=full, aes(x=TOF, y=EXT, color=worm)) + geom_point(alpha=0.5) +
  labs(x="Body length (µm)", y="Optical density") +
  xlim(0, 1500) + ylim(0, 1500) +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        legend.position="none",
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_bubblesROC/bubbleworm.tiff", height=4, width=6, units="in", dpi=300)

#B - ROC curve for difference between worms and bubbles

load(rocplot, file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_bubblesROC/rocplot.RData")

grobroc <- ggplotGrob(ggplot(rocplot) + aes(x=fpr, y=tpr) + geom_line() +
                        xlim(0,0.0075) + ylim(0.85, 1) +
                        labs(x="FPR", y="TPR")+
                        theme_bw() +
                        theme(axis.text.x = element_text(size=10, color="black"),
                              axis.text.y = element_text(size=10, color="black"),
                              axis.title.x = element_text(size=10, face="bold", color="black"),
                              axis.title.y = element_text(size=10, face="bold", color="black"),
                              legend.position="none")
)

ggplot(rocplot) + aes(x=fpr, y=tpr) + geom_line() +
  labs(x="False positive rate", y="True positive rate") +
  annotation_custom(grob=grobroc, xmin=0.25, xmax=1, ymin=0.1, ymax=0.9) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        legend.position="none",
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_bubblesROC/ROC.tiff", height=5, width=5, units="in", dpi=300)

####----------------------------------------------------#####
## Figure S2, incompatibility of QX1430

load("~/Dropbox/HTA + new RIAIL paper/Tables/procincompat.RData")

ggplot(proc) + aes(x=strain, y=m.emb) + 
  geom_bar(position="dodge", stat="identity", fill="darkgray") + 
  geom_errorbar(aes(ymax=m.emb+sd.emb, ymin=m.emb-sd.emb), width=0.35, size=1) +
  labs(x="", y="Mean embryonic lethality") +
  scale_y_continuous(breaks=c(0:5)/10) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, color="black", face="bold", angle=45, vjust=0.65),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_QX1430Incomp/incomp.tiff", height=5, width=6, units="in", dpi=300)

####----------------------------------------------------#####
## Figure S3A
#Plot for allele frequency figure

load("~/Dropbox/HTA + new RIAIL paper/Figures/Figure_AlleleFreq/GenotypeTable_afterImputation.RData")
load("~/Dropbox/HTA + new RIAIL paper/Figures/Figure_AlleleFreq/markers244.Rda")

#convert chromosomes to roman numerals
gt2$chr <- factor(gt2$chr, labels = c(as.character(as.roman(c(1:5, 10)))))

#convert positions to WS.235 positions
gt2$pos <- as.numeric(markers$WS244.pos[match(x=rownames(gt2), markers$SNP)])

ggplot(data=gt2) + aes(x=pos/1000000, y=BB) + geom_line(size=0.5) +
  geom_hline(yintercept=0.5, color="red", linetype=2) +
  facet_grid(.~chr, space="free_x", scale="free_x") +
  xlab("Genomic position (Mb)") + ylab("Frequency of the Bristol allele") + theme_bw() +
  theme(axis.text.x = element_text(size=0, face="bold", color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_AlleleFreq/AF_pts.tiff", height=3, width=6, units="in", dpi=300)

####----------------------------------------------------#####
## Figure S3B
## All new RIAIL genotypes

load(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_newRIAILgenotypes/geno_long.RData")

ggplot(mgen)+
  aes(x = ind, y= strain, color = ifelse(geno=="0","CB4856",
                                         ifelse(geno=="1","N2","2")), 
      fill = ifelse(geno=="0","CB4856",
                    ifelse(geno=="1","N2","2"))) +
  scale_color_manual(values=c("CB4856"="blue","N2"="orange","2"="gray"), name = "Genotype")+
  scale_fill_manual(values=c("CB4856"="blue","N2"="orange","2"="gray"), name = "Genotype")+
  geom_tile(size=1) +
  labs(y = "RIAIL", x = "Marker") +
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=0, face="bold", color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        legend.position="none")

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_newRIAILgenotypes/genopic.tiff", height=6.5, width=6.5, units="in", dpi=300)

####----------------------------------------------------#####
## Figure S4
#Power calculations and B, variance explained

library(qtlDesign)
library(qtl)
library(pwr)

dd = seq(0.01,100,.01)
#n359 = power.t.test(n=180, delta=dd, sig.level=.05/750)$power
n359 = pwr.t.test(n=300, d=dd, sig.level=.05/730)$power #sig corrected for 1460 SNPs
pvar= prop.var('ri', dd/2,1)
df <- data.frame(pvar, n359)

ggplot(data=df) + aes(x=pvar*100, y=n359) + geom_line(size=0.5) + xlim(0, 15) +
  geom_hline(yintercept=0.8, color="red", linetype=2) + 
  xlab("Percent of phenotypic variance explained") + ylab("Statistical power") + theme_bw() +
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_powercalcs/power.tiff", height=3, width=4.5, units="in", dpi=300)



####----------------------------------------------------#####
#Figure S5

fig2DataB <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(!is.na(var.exp))
fig2DataB$chr <- ifelse(fig2DataB$chr==1, "I",
                        ifelse(fig2DataB$chr==2, "II",
                               ifelse(fig2DataB$chr==3, "III",
                                      ifelse(fig2DataB$chr==4, "IV",
                                             ifelse(fig2DataB$chr==5, "V", "X")))))

tabdf2 <- fig2DataB %>% 
  mutate(condition = str_split_fixed(trait, pattern="\\.", n=2)[,1]) %>% 
  mutate(trait2 = str_split_fixed(trait, pattern="\\.", n=2)[,2]) %>%
  select(-trait) %>% rename(trait = trait2) %>%
  filter(!grepl(trait, pattern="f.")) %>%
  filter(!grepl(trait, pattern="react.")) %>%
  arrange(condition, trait, chr)

ggplot(data=tabdf2) + aes(x=var.exp) + geom_histogram() +
  geom_vline(xintercept=0.04, color="red", linetype=2) +
  labs(x="Variance explained", y="Count") +
  theme_bw() + 
  theme(axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_powercalcs/VE.tiff", height=4.5, width=4.5, units="in", dpi=300)

####----------------------------------------------------#####
#Figure -Correlation supplement

pheno <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/SummarizedProcPhenotypes.csv")

#Remove old RIAIL data
pheno <- pheno %>% filter(id>239)

#Separate into control and pq

ctrl <- pheno %>% filter(drug == "control")
pq <- pheno %>% filter(drug=="paraquat")

#Reduce to just control.resid.a and paraquat.resid
ctrl <- ctrl %>% select(id, drug, resid.a.n:resid.a.iqr.norm.EXT)
pq <- pq %>% select(id:resid.iqr.norm.EXT) 

## A

#Melt and look at correlation
cor.ctrl <- melt(cor(ctrl[3:ncol(ctrl)], use="pairwise.complete.obs", method="spearman"))

red.cor.ctrl <- cor.ctrl %>% separate(Var1, into=c("bla", "trait1"), sep="resid.a.") %>% select(-bla) %>%
  separate(Var2, into=c("bla", "trait2"), sep="resid.a.") %>% select(-bla)

red.cor.ctrl$trait1 <- factor(red.cor.ctrl$trait1, levels=c("n", "q10.TOF", "q10.EXT", "q25.TOF", "q25.EXT", 
                                              "median.TOF", "median.EXT", "mean.TOF", "mean.EXT",
                                              "q75.TOF", "q75.EXT", "q90.TOF", "q90.EXT", 
                                              "iqr.TOF", "iqr.EXT", "var.TOF", "var.EXT",
                                              "q10.norm.EXT","q25.norm.EXT","median.norm.EXT","mean.norm.EXT",
                                              "q75.norm.EXT","q90.norm.EXT","iqr.norm.EXT", "var.norm.EXT"))

red.cor.ctrl$trait2 <- factor(red.cor.ctrl$trait2, levels=c("n", "q10.TOF", "q10.EXT", "q25.TOF", "q25.EXT", 
                                                            "median.TOF", "median.EXT", "mean.TOF", "mean.EXT",
                                                            "q75.TOF", "q75.EXT", "q90.TOF", "q90.EXT", 
                                                            "iqr.TOF", "iqr.EXT", "var.TOF", "var.EXT",
                                                            "q10.norm.EXT","q25.norm.EXT","median.norm.EXT","mean.norm.EXT",
                                                            "q75.norm.EXT","q90.norm.EXT","iqr.norm.EXT", "var.norm.EXT"))
#Plot by geom_tile
  
ggplot(red.cor.ctrl)+
  aes(x = trait1, y = trait2)+
  geom_tile(aes(fill=value^2))+
  scale_fill_gradientn(colours=terrain.colors(6), name=expression(bold(italic("r")^2)))+
  theme_bw() +
  theme(axis.text.x = element_text(size=10, face="bold", color="black", angle = 60, hjust=1),
        axis.text.y = element_text(size=10, face="bold", color="black"),
        panel.background = element_rect(fill = '#E5E8E8')) +
  labs(x=NULL, y= NULL)

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_cor/control_corr.tiff", height=7.5, width=7.5, units="in", dpi=300)

## B

#Melt and look at correlation
cor.pq <- melt(cor(pq[3:ncol(pq)], use="pairwise.complete.obs", method="spearman"))

red.cor.pq <- cor.pq %>% separate(Var1, into=c("bla", "trait1"), sep="resid.") %>% select(-bla) %>%
  separate(Var2, into=c("bla", "trait2"), sep="resid.") %>% select(-bla)

red.cor.pq$trait1 <- factor(red.cor.pq$trait1, levels=c("n", "q10.TOF", "q10.EXT", "q25.TOF", "q25.EXT", 
                                                        "median.TOF", "median.EXT", "mean.TOF", "mean.EXT",
                                                        "q75.TOF", "q75.EXT", "q90.TOF", "q90.EXT", 
                                                        "iqr.TOF", "iqr.EXT", "var.TOF", "var.EXT",
                                                        "q10.norm.EXT","q25.norm.EXT","median.norm.EXT","mean.norm.EXT",
                                                        "q75.norm.EXT","q90.norm.EXT","iqr.norm.EXT", "var.norm.EXT"))

red.cor.pq$trait2 <- factor(red.cor.pq$trait2, levels=c("n", "q10.TOF", "q10.EXT", "q25.TOF", "q25.EXT", 
                                                        "median.TOF", "median.EXT", "mean.TOF", "mean.EXT",
                                                        "q75.TOF", "q75.EXT", "q90.TOF", "q90.EXT", 
                                                        "iqr.TOF", "iqr.EXT", "var.TOF", "var.EXT",
                                                        "q10.norm.EXT","q25.norm.EXT","median.norm.EXT","mean.norm.EXT",
                                                        "q75.norm.EXT","q90.norm.EXT","iqr.norm.EXT", "var.norm.EXT"))

#Plot by geom_tile

ggplot(red.cor.pq)+
  aes(x = trait1, y = trait2)+
  geom_tile(aes(fill=value^2))+
  scale_fill_gradientn(colours=terrain.colors(6), name=expression(bold(italic("r")^2)))+
  theme_bw() +
  theme(axis.text.x = element_text(size=10, face="bold", color="black", angle = 60, hjust=1),
        axis.text.y = element_text(size=10, face="bold", color="black"),
        panel.background = element_rect(fill = '#E5E8E8')) +
  labs(x=NULL, y= NULL)

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_cor/pq_corr.tiff", height=7.5, width=7.5, units="in", dpi=300)

#####------------------------------------------------------------------------#####
## Table 1, all QTL mapping results

tabdf <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(!is.na(var.exp))

tabdf2 <- tabdf %>% 
  mutate(condition = str_split_fixed(trait, pattern="\\.", n=2)[,1]) %>% 
  mutate(trait2 = str_split_fixed(trait, pattern="\\.", n=2)[,2]) %>%
  select(-trait) %>% rename(trait = trait2) %>%
  filter(!grepl(trait, pattern="f.")) %>%
  filter(!grepl(trait, pattern="react.")) %>%
  arrange(condition, trait, chr)

write.table(tabdf2, file="~/Dropbox/HTA + new RIAIL paper/Tables/allQTL.csv", sep=",", quote=F, row.names=F, col.names=T)

## Table, all phenotype data for RIAILs

tab2 <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Data/ProcessedPhenotypes.csv")

write.table(tab2, file="~/Dropbox/HTA + new RIAIL paper/Tables/allpheno.csv", sep=",", quote=F, row.names=F, col.names=T)

## Table, all genotype data for new RIAILs, QX240-QX539

genos <- data.frame(fread("~/Dropbox/AndersenLab/LabFolders/Stefan/GWAS/NIL_generation/N2xCB4856_598RIAILs_gen.csv",
                          header = T))

load("~/Dropbox/AndersenLab/RCode/Linkage mapping/markers244.Rda")

genos2 <- do.call(cbind, lapply(genos[,4:ncol(genos)], function(x){
  temp <- data.frame(gen = x)
  temp <- mutate(temp, num = ifelse(gen =="AA",1,
                                    ifelse(gen=="AB",0,NA)))
  temp <- select(temp, -gen)
}))

genos3 <- data.frame(genos[,1:3], genos2)
colnames(genos3) <- c("id","chr","cM", seq(1,ncol(genos3)-3))

genos5 <- genos3 %>%
  gather(key=strain,value=geno,-id  ,-chr,  -cM) %>% 
  mutate(numst = as.numeric(as.character(strain))) %>%
  filter(numst > 239) %>%
  select(-numst) %>%
  spread(key=strain, value=geno) %>%
  arrange(chr, cM)

genos6 <- data.frame(genos5[,1:3], WS244.pos=as.numeric(markers$WS244.pos[match(genos5$id, markers$SNP)]), genos5[,4:362] )

colnames(genos6) <- c("id", "chr", "cM", "WS244.pos", paste0("QX", seq(from=240, to=598)))

#1 = N2, 2 = CB4856

write.table(genos6, file="~/Dropbox/HTA + new RIAIL paper/Tables/allnewgeno.csv", sep=",", quote=F, row.names=F, col.names=T)

#####------------------------------------------------------------------------#####
## Figure - not in paper. Traits are correlated, so hot spots are a little weird to consider.

tabdf <- read.csv("~/Dropbox/HTA + new RIAIL paper/HTA_Linkage/Mapping/MappingResults.csv") %>% filter(!is.na(var.exp))

tabdf2 <- tabdf %>% 
  mutate(condition = str_split_fixed(trait, pattern="\\.", n=2)[,1]) %>% 
  mutate(trait2 = str_split_fixed(trait, pattern="\\.", n=2)[,2]) %>%
  select(-trait) %>% rename(trait = trait2) %>%
  arrange(condition, trait, chr)

tabdf2$grp <- ifelse(tabdf2$condition == "control", 1, 2)

tabdf2$grp <- factor(tabdf2$grp, labels = c("Control", "Paraquat"))

ggplot(tabdf2) + aes(x=pos/1e6) + 
  geom_bar(binwidth=0.5) + 
  labs(x="Position (Mb)", y="Number of QTL") +
  facet_grid(grp~chr, scales="free_y") + 
  theme_bw() +
  theme(axis.text.x = element_text(size=0, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", color="black"),
        strip.text.x = element_text(size=12, face="bold", color="black"),
        strip.text.y = element_text(size=12, face="bold", color="black"),
        plot.title = element_text(size=12, face="bold"))

ggsave(file="~/Dropbox/HTA + new RIAIL paper/Figures/Figure_allLM/allLM.tiff", height=4, width=7.5, units="in", dpi=300)

