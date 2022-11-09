### PS126 microcosms supplementary graphs ###
## Antonia Ahme, 1.11.2022 ##

#### Load packages ####
require(dplyr)
require(vegan)
require(ggplot2)
require(phyloseq)
require(SRS)
require(zCompositions)
require(propr)
require(BiodiversityR)
require(ggrepel)
require(rstatix)

rm(list=ls())
options(stringsAsFactors = F)
setwd("~/AWI/RProjects/PS126")

#### Create color palettes & themes ####
temp_pal <- c("skyblue2", "goldenrod2", "tomato2")

pool_pal <- c(
  "green4",
  "dodgerblue2", 
  "#6A3D9A", # purple
  "#FDBF6F", # lt orange
  "darkturquoise",
  "#FF7F00", # orange
  "skyblue2",
  "deeppink1",
  "#E31A1C", # red
  "green1",
  "maroon", 
  "gold1",
  "orchid1",
  "steelblue4",
  "yellowgreen")

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

pool.theme <- theme(panel.background = element_blank(),
                   panel.border = element_blank(),
                   panel.grid = element_blank(),
                   axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.y=element_text(size= 30, face = "bold"),
                   axis.ticks.y=element_blank(),
                   plot.title = element_blank(),
                   strip.text = element_text(size = 30, face = "bold"),
                   strip.background = element_blank(), strip.placement = "outside",
                   text = element_text(size = 20, face = "bold"), legend.position = "right")


#### Supplementary Figure S1 (Unpooled vs. pooled) MISSING ####
# Import asv and tax tab
asv <- read.delim('Data/Counts_pp.txt')
asv <- data.frame(asv, row.names = 1)
colnames(asv) <- gsub('\\.', '-', colnames(asv))
tax <- read.delim('Data/Taxonomy_pp.txt')

# Import sample tab
sam <- read.delim('Data/Samples_pp.txt')
sam$sample_ID <- as.factor(sam$sample_ID)
sam$temperature <- as.factor(sam$temperature)
sam <- sam %>%
  mutate(temperature = recode(temperature, "2" = '2°C')) %>%
  mutate(temperature = recode(temperature, "6" = '6°C')) %>%
  mutate(temperature = recode(temperature, "9" = '9°C'))

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
rownames(sam) <- sam$sample_ID

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

# Create subsets of temperatures
ps_2 <- subset_samples(ps, temperature == "2°C")
ps_6 <- subset_samples(ps, temperature == "6°C")
ps_9 <- subset_samples(ps, temperature == "9°C")

# Create mean of replicates for temperatures
merged_2 <- merge_samples(ps_2, "bottle", fun = mean)
merged_6 <- merge_samples(ps_6, "bottle", fun = mean)
merged_9 <- merge_samples(ps_9, "bottle", fun = mean)

# Create a table of all classes and abundances
class <- phyloseq::tax_glom(merged_2, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class_2 <- df2 %>% select(Sample, Class, Abundance)
class_2$Temperature <- "2°C"

class <- phyloseq::tax_glom(merged_6, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class_6 <- df2 %>% select(Sample, Class, Abundance)
class_6$Temperature <- "6°C"

class <- phyloseq::tax_glom(merged_9, "Class")
df <- plot_bar(class, fill = "Class")
df2 <- df$data
class_9 <- df2 %>% select(Sample, Class, Abundance)
class_9$Temperature <- "9°C"

class <- rbind(class_2, class_6, class_9)
class <- class %>% rename(Bottle = Sample)

# Rename classes for prettier plotting
class$Class[class$Abundance < 100] <- "Other"
class <- class %>%
  mutate(Class = recode(Class, "MAST-1" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-2" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-3" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-8" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MAST-12" = 'MAST clades')) %>%
  mutate(Class = recode(Class, "MOCH-2" = "MOCH"))
class$Bottle <- factor(class$Bottle, levels = c("Pooled", "Unpooled2", "Unpooled1"))

##Plotting
pool_plot <- ggplot(class, aes(fill = Class, x = Abundance, y = Bottle)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~ Temperature) +
  pool.theme +
  scale_fill_manual(values = pool_pal)

pool_plot
ggsave("Output/SUPP#1.png", pool_plot, height = 8, width = 20, dpi = 320)

#### Supplementary Figure S2 (RDA plot without D2) ####
### Based on the code of Cora Hörstmann ###
## Prepare data
ASV.clr <- read.table('Data/Asv_tfin.txt')
colnames(ASV.clr) <- gsub('\\.', '-', colnames(ASV.clr))
meta <- read.table('Data/Sam_tfin.txt')
meta <- meta[,-c(10)] # remove D2 values

## Create the RDA plot
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

meta_sig <- c('temperature','chl', 'poc', 'cn', 'si', 'rich', 'even')
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
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab(axis.long[1, "label"]) +
  ylab(axis.long[2, "label"]) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_segment(data=subset(vectors.long, vector %in% c("chl", "poc", "si", "cn", "d2", "rich", "even")),
               aes(x=0, y=0, xend=axis1*3, yend=axis2*3), 
               lineend = 'round', linejoin = 'bevel',
               colour="black", size=0.705, arrow=arrow(length = unit(0.2,"cm"))) +
  geom_text_repel(data=subset(vectors.long, vector %in% c("chl", "poc", "si", "cn", "d2", "rich", "even")), 
                  aes(x=axis1*3, y=axis2*3, label=vector), size=5,
                  colour="black") +
  scale_y_reverse() +
  RDA.theme +
  theme(legend.position=c(.8,.85)) +
  coord_fixed(ratio=1)

RDA_plot_cons
ggsave("Output/SUPP#1.png", RDA_plot_cons, height = 6, width = 5.3, dpi = 320)
