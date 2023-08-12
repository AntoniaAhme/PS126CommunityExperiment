### PS126 Bacteria main graphs ###
## Anabel von Jackowski, 24.03.2023##

#### Load packages & housekeeping ####
require(dplyr)
require(vegan)
require(ggplot2)
require(phyloseq)
require(SRS)
require(zCompositions)
require(BiodiversityR) #wierd popup
require(ggrepel)
require(ggpubr)
require(rstatix)
require(easyCODA)
require(edgeR)
require(microbial) #require edgeR
require(ellipse)

rm(list=ls())
options(stringsAsFactors = F)

#### Create color palettes & themes ####
temp_pal <- c("skyblue2", "goldenrod2", "tomato2")

class_pal <- c("green4","dodgerblue2","#FF7F00","darkturquoise","gold1",
               "maroon","deeppink1","darkblue","#E31A1C","peachpuff1","orchid1","aquamarine","darkorange4",
               "steelblue4")

class_pal_35 <- c("yellow","#b4e0ff","#e70017","#00fdc4","#6926d4", "#c0cf00","#b873ff","#0b7d00","#dfbb00","#270044","#e6ff74",
                  "#002664","#b7ff88","#a2006d","#628800","#85aeff","#ff6b31","#00b094","#c30037","#004226","#ff76ab","#1b3f00",
                  "#8f0043","#ffd56f","#020023","#ffac59","#00241f","#ffd8b8","#2d0013","#807000","#6d0021","#ff9c9b","#1f0400",
                  "#9f4d00","#362000")

class_pal_8 <- c(
  "aquamarine4",
  "greenyellow",
  "darkgreen",
  "red",
  "green2",
  "darkolivegreen4",
  "palegreen",
  "forestgreen")

class_pal_4 <- c(
  "gold1",
  "red",
  "yellow1",
  "goldenrod2")

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
asv1 <- read.delim('Data/Counts_bacteria.txt')
asv1 <- data.frame(asv1, row.names = 1)
colnames(asv1) <- gsub('\\.', '-', colnames(asv1))

Cyanobacteria_asv <- subset(tax, Phylum == "Cyanobacteria")
asv <- asv1[!rownames(asv1) %in% rownames(Cyanobacteria_asv),]

tax1 <- read.delim('Data/Taxonomy_bacteria.txt')
tax <- subset(tax1, Phylum !="Cyanobacteria")

# Import sample tab
sam <- read.delim('Data/Samples.txt')
sam$sample_ID <- as.factor(sam$sample_ID)
sam$temperature <- as.factor(sam$temperature)
sam <- sam %>%
  tidyr::unite(uni_ID, timepoint, temperature, sep = "-", remove = FALSE)

# Create subset of sample tabs
sam_tfin <- subset(sam, timepoint == "tfin")
rownames(sam_tfin) <- sam_tfin$sample_ID

#### Prepare data for plotting ####
# Remove ASVs that don't occur in the dataset
asv <- asv[,colnames(asv) %in% sam$sample_ID]

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

# unique(tax$Division)
unique(tax$Class)

# Subset asv tab based on newly selected taxonomy
asv <- asv[rownames(asv) %in% tax$ASV,]

# Create a tax file for raw read phyloseq object
tax_raw <- tax

## Create a phyloseq object of the raw data for diversity estimates
rownames(tax_raw) <- tax_raw$ASV
tax_raw <- tax_raw[,2:7] # delete ASV column
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
tax <- tax[,2:7] # delete ASV column
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

#### Calculate richness and evenness based on ASVs ####
ps.rich <- microbial::richness(ps_raw,  method = c("Observed", "Evenness"))

# Add Shannon diversity to sample tab
sam$Richness <- ps.rich$Observed
sam$Evenness <- ps.rich$Evenness

# Add shannon diversity to t-fin subset of sample tab
sam_tfin <- subset(sam, timepoint == "tfin")
row.names(sam_tfin) <- NULL

# Export sam and asv tab of tfin for supplementary
write.table(sam_tfin, "./Sam_tfin.txt", quote = FALSE, sep= " ")
write.table(asv_tfin, "./Asv_tfin.txt", quote = FALSE, sep= " ")

#### Figure #3 (bargraphs) ####
## Create mean of replicates for t0 and tfin
merged_t0 <- merge_samples(ps_t0, "temperature", fun = mean)
merged_tfin <- merge_samples(ps_tfin, "temperature", fun = mean)

### Figure #3a bargraph on class level of t0 & tfin
## Create a table of all classes and abundances
class <- phyloseq::tax_glom(merged_t0, "Class") # average temperature
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class_t0 <- df2 %>% select(Sample, Class, Abundance)

# class <- phyloseq::tax_glom(ps_tfin, "Class") # no average
class <- phyloseq::tax_glom(merged_tfin, "Class") # average temperature 
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class_tfin <- df2 %>% select(Sample, Class, Abundance) 

class <- rbind(class_t0, class_tfin)
class <- class %>% dplyr::rename(Temperature = Sample)

# Rename classes for prettier plotting
class$Class[class$Abundance < 100] <- "Other"
class$Temperature <- factor(class$Temperature, levels = c("t0", "2C", "6C", "9C"))

# Get percentages of classes (and add manually pre-processing)
class <- aggregate(Abundance~Temperature+Class, data=class, FUN=sum)
maxvalues <- aggregate(Abundance~Temperature, data=class, FUN=sum)
maxvalues

## Classes of PR2 are spanning several taxonomic levels, rename to group
# class <- rename(class, "Class" = "Class")

class_pal <- c("green4","#6A3D9A","darkturquoise","gold1","black",
               "maroon","deeppink1","darkblue","#E31A1C","peachpuff1","orchid1","aquamarine","darkorange4",
               "steelblue4")

## Plotting with and without legend for manual post-processing
classL_plot <- ggplot(class, aes(fill = Class, x = Temperature, y = Abundance)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = class_pal)
classL_plot
class.leg <- get_legend(class_plot)
class.leg <- as_ggplot(class.leg)
class.leg
ggsave("./Experiment_R-Script/Figures/Figure#3ClassLegend.png", class.leg, height = 8, width = 4, dpi = 320)

class_plot <- ggplot(class, aes(fill = Class, x = Temperature, y = Abundance)) +
  geom_bar(position = "fill", stat = "identity") +
  bar.theme +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent")) +
  scale_fill_manual(values = class_pal)

class_plot
ggsave("./Experiment_R-Script/Figures/Figure#3Classes.png", class_plot, height = 5, width = 5, dpi = 320)

### Figure #3b bargraph of Bacteroidia at tfin
## Create a table of all classes and abundances and percent of total community reads
sub <- subset_taxa(merged_t0, Class == "Bacteroidia")
genus <- phyloseq::tax_glom(sub, "Genus")
df <- plot_bar(genus, fill = "Genus")
df2 <- df$data
genus <- df2 %>% select(Sample, Genus, Abundance)
genus <- genus %>% dplyr::rename(Temperature = Sample)
genus$Temperature <- factor(genus$Temperature, levels = c("9C", "6C", "2C"))
percent <- aggregate(Abundance~Temperature, data=genus, FUN=sum)
percent$Max <- as.numeric(percent$Max)
percent$Percent <- (percent$Abundance/percent$Max)*100
BacteroidiaPercent <- percent
BacteroidiaPercent$Genus <- "Bacteroidia"
Bacteroidia

# Rename species for prettier plotting
genus$Genus[genus$Abundance < 50] <- "Other"
genus$Genus <- gsub('_', ' ', genus$Genus)
Bacteroidia <- genus
# write.csv(Bacteroidia, "./Experiment_R-Script/Bacteroidia_Abundance.csv")

## Plotting with and without legend
BacteroidiaL_plot <- ggplot(Bacteroidia, aes(fill = Genus, x = Temperature, y = Abundance)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#0047d5","black","#8100b4","#337495", "red"))
BacteroidiaL_plot
Bacteroidia.leg <- as_ggplot(Bacteroidia.leg)
Bacteroidia.leg
ggsave("./Experiment_R-Script/Figures/Figure#3Bacteroidia_Legend.png", Bacteroidia.leg, height = 6, width = 4, dpi = 320)

Bacteroidia_plot <- ggplot(Bacteroidia, aes(fill = Genus, x = Abundance/sum(class$Abundance), y = Temperature)) +
  geom_bar(position = "stack", stat = "identity") +
  barHor.theme +
  scale_x_continuous(labels = scales::percent, limits = c(0,.21))+
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent")) +
  ggtitle("Bacteroidia") +
  scale_fill_manual(values = c("#0047d5","black","#8100b4","#337495"))

Bacteroidia_plot
ggsave("./Experiment_R-Script/Figures/Figure#3Bacteroidia.pdf", diatoms_plot, height = 5, width = 7, dpi = 320)

### Figure #3c bargraph of haptophytes at tfin
## Create a table of all classes and abundances
sub <- subset_taxa(merged_tfin, Class == "Gammaproteobacteria")
genus <- phyloseq::tax_glom(sub, "Genus")
df <- plot_bar(genus, fill = "Genus")
df2 <- df$data
genus <- df2 %>% select(Sample, Genus, Abundance)
genus <- genus %>% dplyr::rename(Temperature = Sample)
genus$Temperature <- factor(genus$Temperature, levels = c("9C", "6C", "2C"))
percent <- aggregate(Abundance~Temperature, data=genus, FUN=sum)
percent$Max <- as.numeric(percent$Max)
percent$Percent <- (percent$Abundance/percent$Max)*100
GammaPercent <- percent
GammaPercent$Genus <- "Gammaproteobacteria"

# Rename species for prettier plotting
genus$Genus[genus$Abundance < 1000] <- "Other"
genus$Genus <- gsub('_', ' ', genus$Genus)
Gammaproteobacteria <- genus

class_pal_yellow <- c("#fa824a","#6b5d57","#ffd55a","#d08400","#784719","black","#863e06","#ffd290",
                      "#e29978","#b49d62")

## Plotting with and without legend
# write.csv(Gammaproteobacteria, "./Experiment_R-Script/Class_Abundance.csv")
GammaL_plot <- ggplot(Gammaproteobacteria, aes(fill = Genus, x = Temperature, y = Abundance)) +
  geom_bar(position = "stack", stat = "identity") +
  bar.theme +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = class_pal_yellow)
GammaL_plot
Gamma.leg <- get_legend(GammaL_plot)
Gamma.leg <- as_ggplot(Gamma.leg)
Gamma.leg
ggsave("./Experiment_R-Script/Figures/Figure#3GammaLegend.png", Gamma.leg, height = 4, width = 5, dpi = 320)

Gamma_plot <- ggplot(Gammaproteobacteria, aes(fill = Genus, x = Abundance/sum(class$Abundance), y = Temperature)) +
  geom_bar(position = "stack", stat = "identity") +
  barHor.theme +
  scale_x_continuous(labels = scales::percent, limits = c(0,.21))+
  theme(legend.position = "none",
        plot.background = element_rect(fill = "transparent")) +
  ggtitle("Gammaproteobacteria") +
  scale_fill_manual(values = class_pal_yellow)

Gamma_plot
ggsave("./Experiment_R-Script/Figures/Figure#3Gammaproteobacteria.pdf", Gamma_plot, height = 5, width = 7, dpi = 320)

