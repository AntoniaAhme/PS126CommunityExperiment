### PS126 microcosms main graphs ###
## Antonia Ahme, 06.12.2022 ##

#### Load packages & housekeeping ####
require(dplyr)
require(vegan)
require(ggplot2)
require(phyloseq)
require(SRS)
require(zCompositions)
require(propr)
require(BiodiversityR)
require(ggrepel)
require(ggpubr)
require(rstatix)
require(easyCODA)
require(microbial)
require(ellipse)

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/PS126")

#### Create color palettes & themes ####
temp_pal <- c("skyblue2", "goldenrod2", "tomato2")

class_pal <- c(
  "green4",
  "dodgerblue2", 
  "#6A3D9A", # purple
  "skyblue2",
  "#FF7F00", # orange
  "darkturquoise",
  "gold1",
  "maroon", 
  "deeppink1",
  "darkblue",
  "#E31A1C", # red
  "peachpuff1",
  "orchid1",
  "aquamarine",
  "darkorange4",
  "steelblue4")

diatoms_pal <- c(
  "aquamarine4",
  "greenyellow",
  "darkgreen",
  "red",
  "green2",
  "darkolivegreen4",
  "palegreen",
  "forestgreen")

prym_pal <- c(
  "gold1",
  "red",
  "yellow1",
  "goldenrod2")

micro_pal <- c(
  "tomato4",
  "darksalmon",
  "brown3",
  "sienna2")

# Symbols for base R plotting
temp_pch <- c(16,17,15)

# Create the designs for plotting
RDA.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

bar.theme <- theme(panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title.x=element_blank(),
  axis.text.x=element_text(size= 30, face = "bold"),
  axis.ticks.x=element_blank(),
  axis.title.y=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  plot.title = element_blank(),
  strip.text = element_text(size = 30, face = "bold"),
  strip.background = element_blank(), strip.placement = "outside",
  text = element_text(size = 20, face = "bold"), legend.position = "right")

barHor.theme <- theme(panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid = element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_text(size= 30, face = "bold"),
                   axis.ticks.y=element_blank(),
                   plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                   strip.text = element_text(size = 30, face = "bold"),
                   strip.background = element_blank(), strip.placement = "outside",
                   text = element_text(size = 20, face = "bold"), legend.position = "right")


#### Upload & tidy up data ####
# Import asv and tax tab
asv <- read.delim('Data/Counts.txt')
asv <- data.frame(asv, row.names = 1)
colnames(asv) <- gsub('\\.', '-', colnames(asv))
tax <- read.delim('Data/Taxonomy.txt')

# Import sample tab
sam <- read.delim('Data/Samples.txt')
sam$sample_ID <- as.factor(sam$sample_ID)
sam$temperature <- as.factor(sam$temperature)
sam <- sam %>%
  mutate(temperature = recode(temperature, "0" = 't0')) %>%
  mutate(temperature = recode(temperature, "2" = '2°C')) %>%
  mutate(temperature = recode(temperature, "6" = '6°C')) %>%
  mutate(temperature = recode(temperature, "9" = '9°C'))

# Create subset of sample tabs
sam_tfin <- subset(sam, timepoint == "tfin")

#### Prepare data for plotting ####
# Remove ASVs that don't occur in the dataset
asv <- asv[rowSums(asv)>0,]

# Investigate sequencing depth and remove samples where it is too low
depth <- colSums(asv) # prepare df containing the sequencing depth
plot(depth) # sequencing depth does not have any outliers or very low numbers
rm(depth)

# Check and adjust whether rownames of meta-info file match colnames of ASV counts
all(sam$sample_ID == colnames(asv))

# Remove asvs with a count of less than 10 in replicate sample means (create rep means first)
ASV <- asv
colnames(ASV) <- sam$uni_ID[match(colnames(ASV),sam$sample_ID)]
Names <- unique(names(ASV))
ASV <- sapply(Names, function(x)  rowMeans(ASV[names(ASV) %in% x]))
ASV <- as.data.frame(ASV)
ASV <- ASV %>% filter_at(vars(1:4), any_vars(.>=10)) # number of pooled samples and cutoff-level of 10
ASV.rn <- rownames(ASV)
asv <- asv[rownames(asv) %in% ASV.rn,]
rm(ASV,Names,ASV.rn)

# Remove contamination and/or unwanted groups
unique(tax$Division)
unique(tax$Supergroup)
tax <- tax[tax$ASV %in% rownames(asv),] # create matching taxonomy table for e1

# Remove unwanted divisions and NAs in higher taxonomic ranks
tax <- filter(tax, tax$Division!="Metazoa")
tax <- filter(tax, tax$Division!="Fungi")
tax <- filter(tax, tax$Division!="Cryptophyta:nucl")
tax <- tax[!is.na(tax$Division),]

# Rename taxa
tax <- tax %>%
  mutate(Division = recode(Division, "Stramenopiles_X" = 'Stramenopiles')) %>%
  mutate(Division = recode(Division, "Pseudofungi" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "Filosa-Thecofilosea" = 'Thecofilosea')) %>%
  mutate(Class = recode(Class, "Filosa-Imbricatea" = 'Imbricatea')) %>%
  mutate(Class = recode(Class, "Stramenopiles_XX" = 'Stramenopiles')) %>%
  mutate(Class = recode(Class, "Picozoa_X" = 'Picozoa')) %>%
  mutate(Class = recode(Class, "Telonemia_X" = 'Telonemia'))

unique(tax$Division)
unique(tax$Class)

# Subset asv tab based on newly selected taxonomy
asv <- asv[rownames(asv) %in% tax$ASV,]

# Create a tax file for raw read phyloseq object
tax_raw <- tax

## Create a phyloseq object of the raw data for diversity estimates
rownames(tax_raw) <- tax_raw$ASV
tax_raw <- tax_raw[,2:9] # delete ASV column
tax_raw <- as.matrix(tax_raw)
rownames(sam) <- sam$sample_ID
OTU = otu_table(asv, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_raw)
sam$uni_ID <- as.factor(sam$uni_ID)
samples = sample_data(sam)

ps_raw <- phyloseq(OTU, TAX, samples)

## Scaling with ranked subsampling (srs)
# All samples will be scaled to sample with lowest sequencing depth
depth.min <- min(colSums(asv))
asv.srs <- SRS(asv, depth.min) # running the SRS function
rownames(asv.srs) <- rownames(asv)
rm(depth.min)

# Check and adjust that rownames sam tab and colnames of asv tab match
all(rownames(sam) == colnames(asv))

# Prepare tax file for phyloseq object
rownames(tax) <- tax$ASV
tax <- tax[,2:9] # delete ASV column
tax <- as.matrix(tax)

# Name elements for phyloseq object with scaled data
OTU = otu_table(asv.srs, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax)
SAM = sample_data(sam)

# Create phyloseq object
ps <- phyloseq(OTU, TAX, SAM)

# Subset for tfin
ps_t0 <- subset_samples(ps, timepoint=="t0")
ps_tfin <- subset_samples(ps, timepoint!="t0")

## Normalization using CLR transformation
# Remove zeroes
czm <- cmultRepl(t(asv.srs),  label=0, method="CZM") 

# Clr transformation
rho <- propr(czm, metric = 'rho', ivar = 'clr', symmetrize = TRUE, p=0)

# Clr data
clr <- rho@logratio

# Transpose the normalised and filtered asv tab & create dataframe for r
asv.clr <- t(clr)
asv.clr <- as.data.frame(asv.clr)

# Create clr transformed tabs of tfin
asv_tfin <- asv.clr[,colnames(asv.clr) %in% sam_tfin$sample_ID]
tax <- as.data.frame(tax)
tax_tfin <- tax[tax$ASV %in% rownames(asv_tfin),]

#### Calculate richness and evenness based on ASVs ####
ps.rich <- microbial::richness(ps_raw,  method = c("Observed", "Evenness"))

# Add Shannon diversity to sample tab
sam$Richness <- ps.rich$Observed
sam$Evenness <- ps.rich$Evenness

# Add shannon diversity to t-fin subset of sample tab
sam_tfin <- subset(sam, timepoint == "tfin")
row.names(sam_tfin) <- NULL

# Export sam and asv tab of tfin for supplementary
write.table(sam_tfin, "Data/Sam_tfin.txt", quote = FALSE, sep= " ")
write.table(asv_tfin, "Data/Asv_tfin.txt", quote = FALSE, sep= " ")

#### Figure #3 (bargraphs) ####
## Create mean of replicates for t0 and tfin
merged_t0 <- merge_samples(ps_t0, "temperature", fun = mean)
merged_tfin <- merge_samples(ps_tfin, "temperature", fun = mean)

### Figure #3a bargraph on class level of t0 & tfin
## Create a table of all classes and abundances
class <- phyloseq::tax_glom(merged_t0, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class_t0 <- df2 %>% select(Sample, Class, Abundance)

class <- phyloseq::tax_glom(merged_tfin, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class_tfin <- df2 %>% select(Sample, Class, Abundance) 

class <- rbind(class_t0, class_tfin)
class <- class %>% rename(Temperature = Sample)

# Rename classes for prettier plotting
class$Class[class$Abundance < 100] <- "Other"
class <- class %>%
  mutate(Class = recode(Class, "MAST-1" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-2" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-3" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-8" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MOCH-2" = "MOCH")) %>%
  mutate(Class = recode(Class, "Prymnesiophyceae" = "Haptophyta"))
class$Temperature <- factor(class$Temperature, levels = c("t0", "2°C", "6°C", "9°C"))

# Get percentages of classes (and add manually pre-processing)
class <- aggregate(Abundance~Temperature+Class, data=class, FUN=sum)
maxvalues <- aggregate(Abundance~Temperature, data=class, FUN=sum)
maxvalues
#Temperature Abundance
#          t0     88659
#         2°C     88658
#         6°C     88641
#         9°C     88659

## Classes of PR2 are spanning several taxonomic levels, rename to group
class <- rename(class, "Group" = "Class")

## Plotting with and without legend for manual post-processing
classL_plot <- ggplot(class, aes(fill = Group, x = Temperature, y = Abundance)) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = class_pal)
class.leg <- get_legend(classL_plot)
class.leg <- as_ggplot(class.leg)
class.leg
ggsave("Output/Figure#3ClassLegend.png", class.leg, height = 8, width = 4, dpi = 320)

class_plot <- ggplot(class, aes(fill = Group, x = Temperature, y = Abundance)) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  theme(legend.position = "none") +
  scale_fill_manual(values = class_pal)

class_plot
ggsave("Output/Figure#3Classes.png", class_plot, height = 5, width = 5, dpi = 320)

### Figure #3b bargraph of bacillariophytes at tfin
## Create a table of all classes and abundances and percent of total community reads
sub <- subset_taxa(merged_tfin, Class == "Bacillariophyta")
species <- phyloseq::tax_glom(sub, "Species")
df <- plot_bar(species, fill = "Species")
df2 <- df$data
species <- df2 %>% select(Sample, Species, Abundance)
species <- species %>% rename(Temperature = Sample)
species$Temperature <- factor(species$Temperature, levels = c("9°C", "6°C", "2°C"))
percent <- aggregate(Abundance~Temperature, data=species, FUN=sum)
percent <- percent %>%
  mutate(Max = case_when(
    Temperature == '2°C' ~ "88658",
    Temperature == '6°C' ~ "88641",
    Temperature == '9°C' ~ "88659"))
percent$Max <- as.numeric(percent$Max)
percent$Percent <- (percent$Abundance/percent$Max)*100
DiaPercent <- percent
DiaPercent$Species <- "Bacillariophyta"
DiaPercent
species <- species %>%
  mutate(Percent = case_when(
    Temperature == '2°C' ~ "11.6 %",
    Temperature == '6°C' ~ "15.6 %",
    Temperature == '9°C' ~ "58.8 %"))

# Rename species for prettier plotting
species$Species[species$Abundance < 400] <- "Other"
species <- species %>%
  mutate(Species = recode(Species, "Chaetoceros_diadema_1" = 'Chaetoceros_diadema'))
species$Species <- gsub('_', ' ', species$Species)
diatoms <- species

## Plotting with and without legend
diatomsL_plot <- ggplot(diatoms, aes(fill = Species, x = Temperature, y = Abundance)) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  scale_fill_manual(values = diatoms_pal)
dia.leg <- get_legend(diatomsL_plot)
dia.leg <- as_ggplot(dia.leg)
dia.leg
ggsave("Output/Figure#3DiaLegend.png", dia.leg, height = 6, width = 4, dpi = 320)

diatoms_plot <- ggplot(diatoms, aes(fill = Species, x = Abundance, y = Temperature)) +
  geom_bar(position = "stack", stat = "identity") +
  barHor.theme +
  theme(legend.position = "none") +
  ggtitle("Bacillariophyta") +
  scale_fill_manual(values = diatoms_pal)

diatoms_plot
ggsave("Output/Figure#3Diatoms.png", diatoms_plot, height = 5, width = 7, dpi = 320)

### Figure #3c bargraph of haptophytes at tfin
## Create a table of all classes and abundances
sub <- subset_taxa(merged_tfin, Class == "Prymnesiophyceae")
species <- phyloseq::tax_glom(sub, "Species")
df <- plot_bar(species, fill = "Species")
df2 <- df$data
species <- df2 %>% select(Sample, Species, Abundance)
species <- species %>% rename(Temperature = Sample)
species$Temperature <- factor(species$Temperature, levels = c("9°C", "6°C", "2°C"))
percent <- aggregate(Abundance~Temperature, data=species, FUN=sum)
percent <- percent %>%
  mutate(Max = case_when(
    Temperature == '2°C' ~ "88658",
    Temperature == '6°C' ~ "88641",
    Temperature == '9°C' ~ "88659"))
percent$Max <- as.numeric(percent$Max)
percent$Percent <- (percent$Abundance/percent$Max)*100
PrymPercent <- percent
PrymPercent$Species <- "Prymnesiophyceae"

# Rename species for prettier plotting
species <- species %>%
  mutate(Species = recode(Species, "Prymnesiophyceae_Clade_E_XX_sp." = 'Other'))
species$Species <- gsub('_', ' ', species$Species)
prym <- species

## Plotting with and without legend
prymL_plot <- ggplot(prym, aes(fill = Species, x = Temperature, y = Abundance)) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = prym_pal)
prym.leg <- get_legend(prymL_plot)
prym.leg <- as_ggplot(prym.leg)
prym.leg
ggsave("Output/Figure#3PrymLegend.png", prym.leg, height = 3, width = 4, dpi = 320)

prym_plot <- ggplot(prym, aes(fill = Species, x = Abundance, y = Temperature)) +
  geom_bar(position = "stack", stat = "identity") +
  barHor.theme +
  theme(legend.position = "none") +
  ggtitle("Haptophyta") +
  scale_fill_manual(values = prym_pal)

prym_plot
ggsave("Output/Figure#3Prymnesiophytes.png", prym_plot, height = 5, width = 7, dpi = 320)

## Create and export a table of percent of full community to add manually to graphs
Percent <- rbind(DiaPercent, PrymPercent)
require(writexl)
write_xlsx(Percent, "~/AWI/RProjects/PS126/Data/Percent.xlsx")

#### Figure #4 (CA plots) ####
## Based on CODA script of Michael Greenacre 2022

### SIZE 
## Prepare data
size_all <- read.table("Data/Size.txt", header=FALSE)
size_fin <- size_all[, -c(2,3,4)]
group <- size_fin
group <- t(group)
group <- as.data.frame(group)
colnames(group) <- group[1,]
group <- group[-c(1),]
colnames(group)[1] ="Sample"
group[, 5:8] <- sapply(group[, 5:8], as.numeric)

## Prepare and sort the metadata
temp <- as.numeric(factor(group[,3]))
temps <- unique(group[,3])
temps <- sort(temps) 

## Group data with zeroes
group0 <- group[,5:8]
## Number of zeros and percentage of total
sum(group0==0)
# No zeroes in grouped dataset, so no transformation necessary
group <- group0

## Normalizee data
group.n  <- group/rowSums(group)

## Average percentage of size classes
round(100*colMeans(group.n),2)
#uncategorized  Nanoplankton Microplankton  Picoplankton 
#16.70         24.22         55.40          3.68 

## Check if we need weighing
group.n.unw <- CLR(group.n, weight=FALSE)$LR
group.n.unw.var <- apply(group.n.unw, 2, var)
group.n.cm <- colMeans(group.n)
par(mar=c(4.2,4,1,1), font.lab=2, las=1, mfrow=c(1,1))
plot(group.n.cm, group.n.unw.var, log="xy", type="n", 
     xlab="Average compositional values (log-scale)", ylab="Variance of CLR (log-scale)")
text(group.n.cm, group.n.unw.var, labels=colnames(group), col="red", font=4, cex=0.8)
# no weighing necessary

## Show (and sort) the contributions of the groups to the variance
sort(100*group.n.unw.var/sum(group.n.unw.var), decreasing=TRUE)
#Picoplankton  Nanoplankton Microplankton uncategorized 
#56.9116040    38.8738089     3.5291726     0.6854145 

# Boxcox transformation
tgroup <- t(group.n)
tgroup.pro <- CLOSE(tgroup)

tgroup.lra <- LRA(tgroup.pro, weight=FALSE)
tgroup.lra.rpc <- tgroup.lra$rowpcoord 

BoxCox <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- tgroup.pro^alpha
  foo.ca <- CA(foo)
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox[k] <- protest(tgroup.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

## Dimension reduction with CA of square-root transformed compositions
tgroup.ca <- CA(CLOSE(tgroup.pro^0.79))
tgroup.ca$sv <- tgroup.ca$sv/0.79
round(100*tgroup.ca$sv^2/sum(tgroup.ca$sv^2),3)
# percentage contributions of dimensions
#91.406  7.771  0.823

# Prepare data for plotting
tgroup.ca.cpc <- -tgroup.ca$colpcoord
tgroup.ca.rcc <- -tgroup.ca$rowcoord * sqrt(tgroup.ca$rowmass)
tgroup.ca.ctr <- (tgroup.ca.rcc[,1]^2 > 1/nrow(tgroup)) | (tgroup.ca.rcc[,2]^2 > 1/nrow(tgroup)) 
sum(tgroup.ca.ctr)
# Uncategorized is not adding to the solution enough and thus is discarded

## CA biplot of sizes
png(file="Output/Figure#4a.png",width=7,height=5.5,units="in",res=144)
rescale <- 0.4
dim <- c(1,2)
perc.hor <- 91.406; perc.ver <- 7.771
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * rbind(tgroup.ca.cpc, rescale*tgroup.ca.rcc), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
arrows(0, 0, 0.95 * rescale * tgroup.ca.rcc[tgroup.ca.ctr, 1], 
       0.95 * rescale * tgroup.ca.rcc[tgroup.ca.ctr, 2], 
       length = 0.1, angle = 10, col = "black", lwd = 2)
points(tgroup.ca.cpc, pch = temp_pch[temp], col = temp_pal[temp], font = 3, cex = 2)
text(rescale * tgroup.ca.rcc[tgroup.ca.ctr,], labels = rownames(tgroup)[tgroup.ca.ctr], col = "black", 
     cex = 1, font = 2)
set.seed(123)
CIplot_biv(tgroup.ca.cpc[,1], tgroup.ca.cpc[,2], group=temp, groupcols=temp_pal, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           shownames=FALSE)
set.seed(123)
CIplot_biv(tgroup.ca.cpc[,1], tgroup.ca.cpc[,2], group=temp, groupcols=temp_pal, 
           add=TRUE, shade=FALSE, groupnames=temps, alpha=0.99)
dev.off()

### Trophy
## Prepare data
trophy_all <- read.table("Data/Trophy.txt", header=FALSE)
trophy_fin <- trophy_all[, -c(2,3,4)]
group <- trophy_fin
group <- t(group)
group <- as.data.frame(group)
colnames(group) <- group[1,]
group <- group[-c(1),]
colnames(group)[1] ="Sample"
group[, 5:9] <- sapply(group[, 5:9], as.numeric)

## Prepare and sort the metadata
temp <- as.numeric(factor(group[,3]))
temps <- unique(group[,3])
temps <- sort(temps) 

## Group data with zeroes
group0 <- group[,5:9]
## Number of zeros and percentage of total
sum(group0==0)
# No zeroes in grouped dataset, so no transformation necessary
group <- group0

## Normalizee data
group.n  <- group/rowSums(group)

## Average percentage of size classes
round(100*colMeans(group.n),2)
#hetero         photo     parasitic      uncategorized    mixo 
#58.58         33.48          0.79          5.42          1.73  

## Check if we need weighing
group.n.unw <- CLR(group.n, weight=FALSE)$LR
group.n.unw.var <- apply(group.n.unw, 2, var)
group.n.cm <- colMeans(group.n)
par(mar=c(4.2,4,1,1), font.lab=2, las=1, mfrow=c(1,1))
plot(group.n.cm, group.n.unw.var, log="xy", type="n", 
     xlab="Average compositional values (log-scale)", ylab="Variance of CLR (log-scale)")
text(group.n.cm, group.n.unw.var, labels=colnames(group), col="red", font=4, cex=0.8)
# no weighing necessary

## Show (and sort) the contributions of the groups to the variance
sort(100*group.n.unw.var/sum(group.n.unw.var), decreasing=TRUE)
#hetero       parasitic          mixo         photo      uncategorized 
#33.84472      28.43431      20.94184      14.55168       2.22745 

# Boxcox transformation
tgroup <- t(group.n)
tgroup.pro <- CLOSE(tgroup)

tgroup.lra <- LRA(tgroup.pro, weight=FALSE)
tgroup.lra.rpc <- tgroup.lra$rowpcoord 

BoxCox <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- tgroup.pro^alpha
  foo.ca <- CA(foo)
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox[k] <- protest(tgroup.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

## Dimension reduction with CA of square-root transformed compositions
tgroup.ca <- CA(CLOSE(tgroup.pro^0.79))
tgroup.ca$sv <- tgroup.ca$sv/0.79
round(100*tgroup.ca$sv^2/sum(tgroup.ca$sv^2),3)
# percentage contributions of dimensions
#90.862  8.035  1.087  0.016

# Prepare data for plotting
tgroup.ca.cpc <- -tgroup.ca$colpcoord
tgroup.ca.rcc <- -tgroup.ca$rowcoord * sqrt(tgroup.ca$rowmass)
tgroup.ca.ctr <- (tgroup.ca.rcc[,1]^2 > 1/nrow(tgroup)) | (tgroup.ca.rcc[,2]^2 > 1/nrow(tgroup)) 
sum(tgroup.ca.ctr)
# Uncategorized is not adding to the solution enough and thus is discarded

## CA biplot of trophic statuses
png(file="Output/Figure#4b.png",width=7,height=5.5,units="in",res=144)
rescale <- 0.4
dim <- c(1,2)
perc.hor <- 90.862; perc.ver <- 8.035
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * rbind(tgroup.ca.cpc, rescale*tgroup.ca.rcc), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
arrows(0, 0, 0.95 * rescale * tgroup.ca.rcc[tgroup.ca.ctr, 1], 
       0.95 * rescale * tgroup.ca.rcc[tgroup.ca.ctr, 2], 
       length = 0.1, angle = 10, col = "black", lwd = 2)
points(tgroup.ca.cpc, pch = temp_pch[temp], col = temp_pal[temp], font = 3, cex = 2)
text(rescale * tgroup.ca.rcc[tgroup.ca.ctr,], labels = rownames(tgroup)[tgroup.ca.ctr], col = "black", 
     cex = 1, font = 2)
set.seed(123)
CIplot_biv(tgroup.ca.cpc[,1], tgroup.ca.cpc[,2], group=temp, groupcols=temp_pal, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           shownames=FALSE)
set.seed(123)
CIplot_biv(tgroup.ca.cpc[,1], tgroup.ca.cpc[,2], group=temp, groupcols=temp_pal, 
           add=TRUE, shade=FALSE, groupnames=temps, alpha=0.99)
dev.off()

### Thermal niche
## Prepare data
therm_all <- read.table("Data/Therm.txt", header=FALSE)
therm_fin <- therm_all[, -c(2,3,4)]
group <- therm_fin
group <- t(group)
group <- as.data.frame(group)
colnames(group) <- group[1,]
group <- group[-c(1),]
colnames(group)[1] ="Sample"
group[, 5:8] <- sapply(group[, 5:8], as.numeric)

## Prepare and sort the metadata
temp <- as.numeric(factor(group[,3]))
temps <- unique(group[,3])
temps <- sort(temps) 

## Group data with zeroes
group0 <- group[,5:8]
## Number of zeros and percentage of total
sum(group0==0)
# No zeroes in grouped dataset, so no transformation necessary
group <- group0

## Normalizee data
group.n  <- group/rowSums(group)

## Average percentage of size classes
round(100*colMeans(group.n),2)
#Arctic-temperate           Arctic     Cosmopolitan    uncategorized 
#80.93            10.80             5.50             2.78  

## Check if we need weighing
group.n.unw <- CLR(group.n, weight=FALSE)$LR
group.n.unw.var <- apply(group.n.unw, 2, var)
group.n.cm <- colMeans(group.n)
par(mar=c(4.2,4,1,1), font.lab=2, las=1, mfrow=c(1,1))
plot(group.n.cm, group.n.unw.var, log="xy", type="n", 
     xlab="Average compositional values (log-scale)", ylab="Variance of CLR (log-scale)")
text(group.n.cm, group.n.unw.var, labels=colnames(group), col="red", font=4, cex=0.8)
# weighing is necessary, as cosomopolitan with little overall percentage is contributing
# disproportionally much to the variance

## Show (and sort) the contributions of the groups to the variance
sort(100*group.n.unw.var/sum(group.n.unw.var), decreasing=TRUE)
#Cosmopolitan           Arctic Arctic-temperate    uncategorized 
#55.062443        27.257426        13.271006         4.409125 

# symbols for base R plotting
temp_pch <- c(16,17,15)

# Boxcox transformation
tgroup <- t(group.n)
tgroup.pro <- CLOSE(tgroup)

tgroup.lra <- LRA(tgroup.pro, weight=FALSE)
tgroup.lra.rpc <- tgroup.lra$rowpcoord 

BoxCox <- rep(0, 101)
k <- 1
for(alpha in c(seq(1,0.01,-0.01),0.0001)) {
  foo <- tgroup.pro^alpha
  foo.ca <- CA(foo)
  foo.ca.rpc <- foo.ca$rowpcoord 
  BoxCox[k] <- protest(tgroup.lra.rpc, foo.ca.rpc, permutations=0)$t0  
  k <- k+1
}

## Dimension reduction with CA of square-root transformed compositions
tgroup.ca <- CA(CLOSE(tgroup.pro^0.79))
tgroup.ca$sv <- tgroup.ca$sv/0.79
round(100*tgroup.ca$sv^2/sum(tgroup.ca$sv^2),3)
# percentage contributions of dimensions
#89.043  7.320  3.637

# Prepare data for plotting
tgroup.ca.cpc <- -tgroup.ca$colpcoord
tgroup.ca.rcc <- -tgroup.ca$rowcoord * sqrt(tgroup.ca$rowmass)
tgroup.ca.ctr <- (tgroup.ca.rcc[,1]^2 > 1/nrow(tgroup)) | (tgroup.ca.rcc[,2]^2 > 1/nrow(tgroup)) 
sum(tgroup.ca.ctr)
# Uncategorized is not adding to the solution enough and thus is discarded

## CA biplot of thermal niches
png(file="Output/Figure#4c.png",width=7,height=5.5,units="in",res=144)
rescale <- 0.4
dim <- c(1,2)
perc.hor <- 89.043; perc.ver <- 7.320
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8)
plot(1.05 * rbind(tgroup.ca.cpc, rescale*tgroup.ca.rcc), type = "n", asp = 1, 
     xlab = paste("CA dimension ", dim[1], " (", round(perc.hor, 1), "%)", sep = ""), 
     ylab = paste("CA dimension ", dim[2], " (", round(perc.ver, 1), "%)", sep = ""), 
     xaxt = "n", yaxt = "n", main = "")
abline(h = 0, v = 0, col = "gray", lty = 2)
axis(1)
axis(2)
arrows(0, 0, 0.95 * rescale * tgroup.ca.rcc[tgroup.ca.ctr, 1], 
       0.95 * rescale * tgroup.ca.rcc[tgroup.ca.ctr, 2], 
       length = 0.1, angle = 10, col = "black", lwd = 2)
points(tgroup.ca.cpc, pch = temp_pch[temp], col = temp_pal[temp], font = 3, cex = 2)
text(rescale * tgroup.ca.rcc[tgroup.ca.ctr,], labels = rownames(tgroup)[tgroup.ca.ctr], col = "black", 
     cex = 1, font = 2)
set.seed(123)
CIplot_biv(tgroup.ca.cpc[,1], tgroup.ca.cpc[,2], group=temp, groupcols=temp_pal, 
           add=TRUE, shade=TRUE, alpha=0.99, 
           shownames=FALSE)
set.seed(123)
CIplot_biv(tgroup.ca.cpc[,1], tgroup.ca.cpc[,2], group=temp, groupcols=temp_pal, 
           add=TRUE, shade=FALSE, groupnames=temps, alpha=0.99)
dev.off()

#### Figure #6 (RDA plot) ####
### Based on the code of Cora Hörstmann ###
## Prepare data
ASV.clr <- asv_tfin
ASV.clr <- ASV.clr[, -c(9)] # remove 6°C C due to missing D2
meta <- sam_tfin
meta <- meta[-c(9),] # remove 6°C C due to missing D2

## Create RDA plot
# ASV data prep (sorting according to metadata)
ASV.clr.t <- t(ASV.clr)
ASV.clr.t.sort <- as.data.frame(ASV.clr.t)
ASV.clr.t.sort <- cbind(ASV.clr.t.sort, meta$sample_ID)
ASV.clr.t.sort <- with(ASV.clr.t.sort, ASV.clr.t.sort[order(meta$sample_ID),])
ASV.clr.t.sort$`meta$sample_ID` <-NULL
ASV.clr.t.sort <- as.matrix(data.matrix(ASV.clr.t.sort))
meta.wf <- with(meta, meta[order(sample_ID),])
meta.wf$uni_ID <- NULL
meta.wf$replicate <- NULL
meta.wf$sample_ID <- NULL
meta.wf$timepoint <- NULL
meta.wf <-rename(meta.wf, "C:N" = "C.N")

meta_sig <- c('temperature','Chla', 'POC', 'C:N', 'Silicate', 'D2', 'Richness', 'Evenness')
meta.wf.data <- meta.wf[meta_sig]
  
# Perform RDA
ASV.clr.rda <- rda(
    ASV.clr.t.sort ~ .,
    data = meta.wf.data,
    na.action = na.omit)
  
print(ASV.clr.rda)

# Add metadata as constraints
Condit.env <- meta.wf

invisible(hist(residuals(ASV.clr.rda), main = ""))
  
ASV.clr.rda.anova <- anova.cca(ASV.clr.rda)
print(ASV.clr.rda.anova)
  
inertia.rda.tot <- ASV.clr.rda$tot.chi
inertia.rda.tot
inertia.rda.constrained <- ASV.clr.rda$CCA$tot.chi
inertia.rda.constrained
inertia.rda.constrained.prop <- inertia.rda.constrained/inertia.rda.tot
print(inertia.rda.constrained.prop)
  
# Extract the first two RDA axes for plotting
Ordination.model <- ASV.clr.rda
attach(Condit.env)
summary(Ordination.model)
plot <- ordiplot(Ordination.model, choices=c(1,2))
  
sites.long <- sites.long(plot, env.data=Condit.env)
head(sites.long)
  
axis.long <- axis.long(Ordination.model, choices=c(1, 2))
axis.long
  
vectors.envfit <- envfit(plot, env=Condit.env)
vectors.long <- vectorfit.long(vectors.envfit)
vectors.long

names(sites.long)[names(sites.long) == 'temperature'] <- 'Temperature'

RDA_plot_cons <- ggplot(data=sites.long, 
                    aes(x=axis1, y=axis2, color=Temperature)) + 
    geom_point(size=3, aes(shape = Temperature, colour = Temperature)) +
    scale_color_manual(values = temp_pal) +
    scale_fill_manual(values = temp_pal) +
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axis.long[1, "label"]) +
    ylab(axis.long[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_segment(data=subset(vectors.long, vector %in% c('Chla', 'POC', 'C:N', 'Silicate', 'D2', 'Richness', 'Evenness')),
                 aes(x=0, y=0, xend=axis1*3, yend=axis2*3), 
                 lineend = 'round', linejoin = 'bevel',
                 colour="black", size=0.705, arrow=arrow(length = unit(0.2,"cm"))) +
    geom_text_repel(data=subset(vectors.long, vector %in% c('Chla', 'POC', 'C:N', 'Silicate', 'D2', 'Richness', 'Evenness')), 
                    aes(x=axis1*3, y=axis2*3, label=vector), size=5,
                    colour="black") +
#    ggforce::geom_mark_ellipse(aes(fill = Temperature, color = Temperature)) +
    RDA.theme +
    theme(legend.position=c(.8,.85)) +
    coord_fixed(ratio=1)
  
RDA_plot_cons
ggsave("Output/Figure#6.png", RDA_plot_cons, height = 6, width = 5.3, dpi = 320)

#### Statistics of functional parameters ####
## Upload data
require(car)
fun <- read.table("Data/Sam_tfin.txt", header = TRUE)
fun$sample_ID <- NULL
fun$timepoint <- NULL
fun$uni_ID <- NULL
fun$temperature <- as.factor(fun$temperature)
fun <- rename(fun, "CN" = "C.N")

# Check all parameters for normality, symmetry and homogeneity of variance
boxplot(fun$Chla ~ fun$temperature)
leveneTest(fun$Chla, fun$temperature)
#exchange parameter one after the other
#log-transform all data first, as at least one temp per parameter is skewed
#variances are fine

## Perform bonferroni-corrected t-tests
# chlorophyll
fun$Chla <- log(fun$Chla)
chl <- fun %>%
  pairwise_t_test(Chla ~ temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Chla ~ temperature, data = fun)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# POC
fun$POC <- log(fun$POC)
poc <- fun %>%
  pairwise_t_test(POC ~ temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(POC ~ temperature, data = fun)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# CN
fun$CN <- log(fun$CN)
cn <- fun %>%
  pairwise_t_test(CN ~ temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(CN ~ temperature, data = fun)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# si
fun$Silicate <- log(fun$Silicate)
si <- fun %>%
  pairwise_t_test(Silicate ~ temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Silicate ~ temperature, data = fun)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# d2
fun$D2 <- log(fun$D2)
d2 <- fun %>%
  pairwise_t_test(D2 ~ temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(D2 ~ temperature, data = fun)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# rich
fun$Richness <- log(fun$Richness)
rich <- fun %>%
  pairwise_t_test(Richness ~ temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Richness ~ temperature, data = fun)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# even
fun$Evenness <- log(fun$Evenness)
even <- fun %>%
  pairwise_t_test(Evenness ~ temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Evenness ~ temperature, data = fun)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

## Create a dataframe with significances
pairs <- c("2°C-6°C", "2°C-9°C", "6°C-9°C")
sig <- data.frame(pairs)

sig$chl <- chl$p.adj
sig$poc <- poc$p.adj
sig$cn <- cn$p.adj
sig$si <- si$p.adj
sig$d2 <- d2$p.adj
sig$rich <- rich$p.adj
sig$even <- even$p.adj

require(writexl)
write_xlsx(sig, "Data/FunctionsSignificances.xlsx")

#### Statistics of traits ####
## Upload data
require(car)
sista <- read.table("Data/SizeStat.txt", header = TRUE)
sista$Sample <- NULL
sista$Temperature <- as.factor(sista$Temperature)

trosta <- read.table("Data/TrophyStat.txt", header = TRUE)
trosta$Sample <- NULL
trosta$Temperature <- as.factor(trosta$Temperature)

thesta <- read.table("Data/ThermStat.txt", header = TRUE)
thesta$Sample <- NULL
thesta$Temperature <- as.factor(thesta$Temperature)

nine <- subset(trosta, Temperature == "9")
boxplot(nine$Para ~ nine$Replicate)

# Check all parameters for normality and symmetry
boxplot(sista$Pico ~ sista$Temperature)
boxplot(sista$Nano ~ sista$Temperature)
boxplot(sista$Micro ~ sista$Temperature)
boxplot(trosta$Hetero ~ trosta$Temperature)
boxplot(trosta$Photo ~ trosta$Temperature)
boxplot(trosta$Mixo ~ trosta$Temperature)
boxplot(trosta$Para ~ trosta$Temperature)
boxplot(thesta$Arctic ~ thesta$Temperature)
boxplot(thesta$Arctic.temperate ~ thesta$Temperature)
boxplot(thesta$Cosmopolitan ~ thesta$Temperature)
#log-transform all data first, as at least one temp per parameter is skewed

leveneTest(sista$Pico ~ sista$Temperature)
leveneTest(trosta$Hetero ~ trosta$Temperature)
leveneTest(thesta$Arctic ~ thesta$Temperature)
#exchange parameter one after the other -> all have common variances


## Perform bonferroni-corrected t-tests
# Picoplankton
sista$Pico <- log(sista$Pico)
pico <- sista %>%
  pairwise_t_test(Pico ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Pico ~ Temperature, data = sista)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Nanoplankton
sista$Nano <- log(sista$Nano)
nano <- sista %>%
  pairwise_t_test(Nano ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Nano ~ Temperature, data = sista)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Microplankton
sista$Micro <- log(sista$Micro)
micro <- sista %>%
  pairwise_t_test(Micro ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Micro ~ Temperature, data = sista)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Heterotrophs
trosta$Hetero <- log(trosta$Hetero)
hetero <- trosta %>%
  pairwise_t_test(Hetero ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Hetero ~ Temperature, data = trosta)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Phototrophs
trosta$Photo <- log(trosta$Photo)
photo <- trosta %>%
  pairwise_t_test(Photo ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Photo ~ Temperature, data = trosta)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Mixotrophs
trosta$Mixo <- log(trosta$Mixo)
mixo <- trosta %>%
  pairwise_t_test(Mixo ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Mixo ~ Temperature, data = trosta)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Parasites
trosta$Para <- log(trosta$Para)
para <- trosta %>%
  pairwise_t_test(Para ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Para ~ Temperature, data = trosta)
summary(res.aov)
TukeyHSD(res.aov)
# Here, there are significant differences between 9 and 2 degrees

# Arctic species
thesta$Arctic <- log(thesta$Arctic)
arc <- thesta %>%
  pairwise_t_test(Arctic ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Arctic ~ Temperature, data = thesta)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Arctic-temperate species
thesta$Arctic.temperate <- log(thesta$Arctic.temperate)
tem <- thesta %>%
  pairwise_t_test(Arctic.temperate ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Arctic.temperate ~ Temperature, data = thesta)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

# Cosmopolitan species
thesta$Cosmopolitan <- log(thesta$Cosmopolitan)
cos <- thesta %>%
  pairwise_t_test(Cosmopolitan ~ Temperature, p.adjust = "bonferroni")
# Compare to one-way ANOVA with post-hoc TUKEY t-tests
res.aov <- aov(Cosmopolitan ~ Temperature, data = thesta)
summary(res.aov)
TukeyHSD(res.aov)
# Results do not deviate from pairwise a priori t-tests

## Create a dataframe with significances
pairs <- c("2°C-6°C", "2°C-9°C", "6°C-9°C")
sigTr <- data.frame(pairs)

sigTr$Pico <- pico$p.adj
sigTr$Nano <- nano$p.adj
sigTr$Micro <- micro$p.adj
sigTr$Hetero <- hetero$p.adj
sigTr$Photo <- photo$p.adj
sigTr$Mixo <- mixo$p.adj
sigTr$Para <- para$p.adj
sigTr$Arctic <- arc$p.adj
sigTr$Temperate <- tem$p.adj
sigTr$Cosmo <- cos$p.adj

require(writexl)
write_xlsx(sigTr, "Data/TraitsSignificances.xlsx")

