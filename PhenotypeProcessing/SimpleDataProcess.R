library(COPASutils)
library(plyr)
library(dplyr)
library(kernlab)
library(knitr)

# Sources the required outside functions

source("PhenotypeProcessing/SimpleDataProcessFxns.R")

# Determines whether reports should be generated or not (TRUE to generate reports, FALSE to skip report generation)
# Note: Report generation takes ~45 min

generateReports <- FALSE

# Set all experiment directories

directories <- c("Data/20110811_RIAILs0a", "Data/20110812_RIAILs0b", "Data/20110818_RIAILs0c", "Data/20110819_RIAILs0d")

#Create the folders for the reports in eachg experiment's data folder

sapply(directories, function(x){dir.create(file.path(x, "reports"), showWarnings = FALSE)})

# Get the vector of setup files

setupList <- sapply(unlist(sapply(directories,
                                  function(x){lapply(list.files(x, pattern="\\.txt", recursive=TRUE, full.names=TRUE),
                                                     function(y){if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("setup",y)){return(y)}})})
), as.character)

# Get the vector of score files
scoreList <- sapply(unlist(sapply(directories,
                                  function(x){lapply(list.files(x, pattern="\\.txt", recursive=TRUE, full.names=TRUE),
                                                     function(y){if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("score",y)){return(y)}})})
), as.character)

# Read in all setup plates

sortData <- do.call(rbind, lapply(setupList, function(x){data.frame(cbind(info(x,2), procSetup(x, 60, 2000), step="setup"))}))

# Generate reports, if required

if(generateReports){
    sortData %>% group_by(assay, plate) %>% do(data.frame(setupReport(.)))
}

# Read in all score plates

rawScoreData <- data.frame(do.call(rbind, lapply(scoreList, function(x){data.frame(cbind(info(x,2), readPlate_worms(x, 60, 2000), step="score"))}))) %>% filter(drug != "missing") %>% data.frame()

#save(rawScoreData, file="Data/RawScoreData.Rda")

# Get the strain data

strainFiles <- do.call(rbind, lapply(lapply(directories,
                                            function(x){
                                                lapply(list.files(x,
                                                                  pattern="\\.Rds",
                                                                  recursive=TRUE,
                                                                  full.names=TRUE),
                                                       function(y){
                                                           if(!grepl("IncompleteData",y) & !grepl("UnstitchedData",y) & !grepl("Archive",y) & grepl("strains",y)){
                                                               return(y)
                                                           }
                                                       })
                                            }), function(x){
                                                y=unlist(x)
                                                data.frame(cbind(y[1], y[2]))
                                            }))

strainsData <- do.call(rbind, apply(strainFiles, 1, function(x){
    load(x[1])
    load(x[2])
    ctrls <- data.frame(assay=info(x[1])$assay, plate=1:length(ctrl.strains), strains=I(ctrl.strains))
    pq <- data.frame(assay=info(x[2])$assay, plate=((length(ctrl.strains)+1):(length(ctrl.strains)+length(pq.strains))), strains=I(pq.strains))
    return(data.frame(rbind(ctrls, pq)))
}))
colnames(strainsData) <- c("assay", "plate", "strains")

# Summarize the score data

summarizedScoreData <- rawScoreData %>%
    group_by(date, experiment, round, assay, plate, drug) %>%
    do({summarizePlate_worms(., strains=eval(strainsData[as.character(strainsData$assay)==as.character(.$assay[1]) & as.numeric(strainsData$plate)==as.numeric(as.character(.$plate[1])), "strains"])[[1]], quantiles=TRUE)})

# Join the number sorted and calculate norm.n, then remove the number sorted

completeData <- left_join(summarizedScoreData, select(sortData, assay, plate, row, col, n.sorted.setup = n.sorted)) %>% mutate(norm.n=n/n.sorted.setup) %>% select(-contains("sort"))

# NA out wash wells

completeData[is.na(completeData$strain) | is.na(completeData$norm.n), which(colnames(completeData)=="n"):ncol(completeData)] <- NA

# Generate the score reports

if(generateReports){
    completeData %>% group_by(assay, plate) %>% do(data.frame(scoreReport(.)))
    unlink("figure", recursive=TRUE)
}

# Get a data frame of all of the control files (which plates match up with which control plates)

controls <- data.frame()

for(i in unique(completeData$assay)){
    data <- filter(completeData, assay==i)
    controlPlates <- unique(as.numeric(as.character(filter(data, drug=="control")$plate)))
    pqPlates <- unique(as.numeric(as.character(filter(data, drug=="paraquat")$plate)))
    row <- data.frame(matrix(nrow=1, ncol=3))
    row[[1,1]] <- i
    row[[1,2]] <- parse(text=as.character(call("c", controlPlates))[2])
    row[[1,3]] <- parse(text=as.character(call("c", pqPlates))[2])
    controls <- rbind(controls, row)
}

colnames(controls) <- c("assay", "control", "plates")

# Calculate residuals and reaction norms for all traits

finalData <- completeData %>% group_by(drug) %>% do(regress(., completeData, controls)) %>% arrange(assay) %>% data.frame()
#finalData[finalData$drug=="control",which(colnames(finalData)=="resid.n"):ncol(finalData)] <- NA

ggplot(finalData) + aes(x=n) + geom_histogram() + facet_grid(drug~.)

#Need to prune outlier wells from data

#Remove outlier brood size wells
red.final <- finalData %>% 
  filter(n > 5) %>%
  filter(n<330) %>%
  filter(!is.na(strain))

ggplot(red.final) + aes(x=n) + geom_histogram() + facet_grid(drug~.)

#Reduce data to resid and resid.a. Also, eliminate all f. traits.

red.final2 <- red.final %>% select(assay:strain, resid.n:resid.iqr.yellow, resid.a.n:resid.a.iqr.yellow) %>% data.frame()

write.csv(red.final2, file="Data/ProcessedPhenotypes.csv", row.names=FALSE)
