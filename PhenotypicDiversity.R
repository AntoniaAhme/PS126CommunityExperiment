#### PS126 calculating phenotypic diversity from flow cytometric counts ####
### Antonia Uthoff, 04.10.2022
## On basis of http://rprops.github.io/PhenoFlow/

#### Housekeeping ####
rm(list=ls())

library(dplyr)
library(tidyr)
library(rlang)
library(flowCore)
library(lattice)
library(flowViz) 
library(ggplot2)
library(flowAI)
library(flowFP)
library(flowFDA)
library(car)
library(Phenoflow)
library(gam)
library(gridExtra)
library(grid)
library(writexl)

set.seed(777)

#### Upload, tidy up and transform the data ####
path = "C:/Users/authoff/Documents/AWI/RProjects/PS126/Data/FCS"
flowData <- read.flowSet(path = path, pattern=".fcs")

attributes(flowData)

### Extract metadata from sample names and merge
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData),"_"), rbind)))
colnames(metadata) <- c("timepoint", "temp", "rep")
metadata %>%
  unite("sampleNames", c("timepoint", "temp", "rep"), sep = "_", remove = FALSE) %>%
  unite("ID", c("timepoint","temp"), sep = "_", remove = FALSE) %>%
  separate(rep, c("rep","leftover"), sep = -4) -> metadata
metadata$leftover <- NULL
rownames(metadata) <- metadata$sampleNames
metadata$sampleNames <- NULL

# Add it to flowSet and label it
pData <- as(metadata, "AnnotatedDataFrame")
phenoData(flowData) <- pData
varMetadata(phenoData(flowData))[, "labelDescription"] <-
  c("ID", "Timepoint", "Temperature", "Replicate", "Filename")

### First transformation of flowSet for visualisation
flowData_transformed <- transform(flowData,`FL3-H`=asinh(`FL3-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL2-H`=asinh(`FL2-H`),
                                  `FSC-H`=asinh(`FSC-H`), 
                                  `FL4-H`=asinh(`FL4-H`))
param=c("FL3-H", "SSC-H","FL2-H","FSC-H", "FL4-H")

### Remove the bead data
path = "C:/Users/authoff/Documents/AWI/RProjects/PS126/Data"
flowDataBlank <- read.flowSet(path = path, pattern=".fcs")

attributes(flowDataBlank)

flowDataBlank_transformed <- transform(flowDataBlank,`FL3-H`=asinh(`FL3-H`), 
                                       `SSC-H`=asinh(`SSC-H`), 
                                       `FL2-H`=asinh(`FL2-H`),
                                       `FSC-H`=asinh(`FSC-H`), 
                                       `FL4-H`=asinh(`FL4-H`))
param=c("FL3-H", "SSC-H","FL2-H","FSC-H", "FL4-H")

sqrcut1 <- matrix(c(11.25,11.25,12.5,12.5,13,14.8,14.8,13),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL3-H","SSC-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Beads")

###  Gating quality check
# In Blank
xyplot(`FL3-H` ~ `SSC-H`, data=flowDataBlank_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(7,16)),
                   x=list(limits=c(0,17))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

# In real dataset
xyplot(`FL3-H` ~ `SSC-H`, data=flowData_transformed[1], filter=polyGate1,
       scales=list(y=list(limits=c(7,16)),
                   x=list(limits=c(0,17))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

### Create a subset with and without the bead data based on the blank
flowData_transformed_Beads <- filter(flowData_transformed, polyGate1)
Subset(flowData_transformed, flowData_transformed_Beads)
flowData_1 <- split(flowData_transformed, flowData_transformed_Beads)

###  Gating quality check
xyplot(`FL3-H` ~ `SSC-H`, data=flowData_1$`Beads-`[1], filter=polyGate1,
       scales=list(y=list(limits=c(7,16)),
                   x=list(limits=c(0,17))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

# Tidy up
flowData_final <- flowData_1$`Beads-`
rm(flowData_1, flowData_transformed, flowData_transformed_Beads, flowDataBlank, flowDataBlank_transformed, polyGate1, sqrcut1)

### Normalise to size (SSC-H) counts
summary <- fsApply(x = flowData_final, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
maxval <- max(summary[,"SSC-H"]) 
mytrans <- function(x) x/maxval
flowData_final <- transform(flowData_final,`FL3-H`=mytrans(`FL3-H`),
                            `FL2-H`=mytrans(`FL2-H`),
                            `SSC-H`=mytrans(`SSC-H`),
                            `FSC-H`=mytrans(`FSC-H`),
                            `FL4-H`=mytrans(`FL4-H`))

### Randomly resample to lowest sample size
flowData_final <- FCS_resample(flowData_final)

#### Calculate, export and visualise phenotypic diversity ####
### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_final, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)


### Calculate Diversity and export as xlsx
Diversity.fbasis <- Diversity(fbasis, d=3, plot=FALSE, R=999)

write_xlsx(Diversity.fbasis,"C:/Users/authoff/Documents/AWI/RProjects/PS126/Data/PhenotypicDiversity.xlsx")

### Plot D2 of phenotypic finterprint
levels(metadata$temp) <- c("0", "2", "6", "9")
tempPalette <- c("0" = "grey", "2" = "skyblue2", "6" = "goldenrod2", "9" = "tomato2")

p1 <- ggplot(data = Diversity.fbasis, aes(x = metadata$timepoint, y = D2), group = metadata$ID) +
  geom_point(size = 10, alpha = 0.8, aes(color = metadata$temp)) +
  theme_classic() +
#  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), size=1.1, width=.1,
#                position=position_dodge(0.05)) +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 15)) +
  labs(y = "Phenotypic diversity (D2)", x = "Timepoint") +
  scale_color_manual(name = "Temperature [°C]", values = tempPalette)
p1

p2 <- qplot(metadata$timepoint, D2, data = Diversity.fbasis, 
      geom = "boxplot", fill = metadata$temp) +
      theme_classic() +
      theme(axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 15)) +
      labs(y = "Phenotypic diversity (D2)", x = "Timepoint") +
      scale_fill_manual(name = "Temperature [°C]", values = tempPalette)
p2
