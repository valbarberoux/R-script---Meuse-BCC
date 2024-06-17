#Path####
path<-"/Users/valentinbarberoux/Desktop/R/Donnees_GenoToul/Seq_poolés/Total"
setwd(dir=path)
list.files(path)

library(dada2)
library(phyloseq) 
library (magrittr)
library(ggplot2)
library(vegan)
library(effects)
library(plyr)
library(devtools)
library(ggplot2)
library(htmlwidgets)
library(plotly)
library(dplyr)
library(breakaway)
library(ape)
library(data.table)
library(RColorBrewer)
library(reshape)
library(reshape2)
library(BiocParallel)
library(tidyverse)
library(gtable)
library(grid)
library(gridExtra)
library(ggrepel)
library(adespatial)
library(metagenomeSeq)
library(ALDEx2)
library(DESeq2)
library(ggvegan)
library(DivNet)
library(metagMisc)

# DADA2 pipeline####

fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

length(fnFs)
length(fnRs)

# Place filtered files in filtered/ subdirectory


filtFs <- file.path(path, "filtered", paste0(sample.names, "_R1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

plotQualityProfile(fnFs[1:96])

plotQualityProfile(fnRs[1:96])

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,210), ## trunc the length of the sequences when the sequencing quality start to have a low quality
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft= c(19,20)) # That is for Mac, on Windows set multithread=FALSE
head(out)



errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample and the third one for instance
head(mergers[[1]])
head(mergers[[3]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths 
table(nchar(getSequences(seqtab)))


seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 251:256]
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getX<-getSequences(seqtab.nochim)[1:5]

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) 

# Assign taxonomy using silva database
taxa <- assignTaxonomy(seqtab.nochim, "/Users/valentinbarberoux/Desktop/R/Last silva database/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Export OTU table from seqtabnochim####
library("writexl")
seqtab.nochim<-t(seqtab.nochim) # necessary to tranpose to be able to open it in excel, otherwise not possible
df<-data.frame(seqtab.nochim)
write_xlsx(df,"OTU_tableseqtabnochim") 

# sequences details of each ASV are exported too in this file
seqtab.nochim<-t(seqtab.nochim)
df<-data.frame(seqtab.nochim)
write_xlsx(df,"OTU_table")

# Export TAX table from "taxa"####
df<-data.frame(taxa)
write_xlsx(df,"TAX_table")




# Random rarefaction#####

# Necessary before performing all alpha diversity analyses, but not for the beta diversity analyses as a log trnasfomration is used for that
A <- seqtab.nochim
C <- colSums(A)
min(colSums(A))# is set at 10019 even if some samples had less reads, in order not to lose too much diversity information
B <- rrarefy((A), sample = 10019)
# otherwise, if no manual threshold has to be set, the following command can be used
Bbis <- rrarefy((A), sample = min(colSums(A)))


# Export OTU table from seqtabnochim rarefied####

library("writexl")
df<-data.frame(B)
write_xlsx(df,"OTU_table_rarefied") # sequences pas exportées, pour les exporter transposer (voir ci-dessous)

# Creation of the PSEQ object with rarefied OTU table ###

OTU<-(read.delim2("OTU_table_rarefied.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) # xlxs file are always transormed in txt file before being loaded in R
TAX<-as.matrix(read.delim("TAX_table.txt",row.names=1,header=TRUE)) # correspond à mon tax_table que j'ai exporté après DADA2 pipeline
env_data <- sample_data(read.delim2("metadata.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)


OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
PSEQ <- phyloseq(OTU,TAX,env_data)
PSEQ


# Remove unwanted taxa ####
PSEQ_unclassified_excluded0= subset_taxa(PSEQ, !Kingdom=='NA')
PSEQ_unclassified_excluded0
PSEQ_unclassified_excluded1= subset_taxa(PSEQ_unclassified_excluded0, !Kingdom=='Eukaryota_unclassified')
PSEQ_unclassified_excluded1
PSEQ_unclassified_excluded2= subset_taxa(PSEQ_unclassified_excluded1, !Kingdom=='Archaea')
PSEQ_unclassified_excluded2
PSEQ_unclassified_excluded3= subset_taxa(PSEQ_unclassified_excluded2, !Order=='Chloroplast')
PSEQ_unclassified_excluded3
PSEQ_unclassified_excluded4= subset_taxa(PSEQ_unclassified_excluded3, !Family=='Mitochondria')

PSEQ_last_rarefied<-PSEQ_unclassified_excluded4

# export otu_table  et tax table without unwanted taxa
df<-data.frame(otu_table (PSEQ_last_rarefied))
write_xlsx(df,"PSEQ_withoutunwantedtaxa_otu_table")
df<-data.frame(tax_table (PSEQ_last_rarefied))
write_xlsx(df,"PSEQ_withoutunwantedtaxa_tax_table")


# Make subset of samples####

Temporal<- subset_samples(PSEQ_last_rarefied,Typeofstudy=='Temporal')
Temporal_SF<- subset_samples(Temporal,Fraction=='SF')
Temporal_LF<- subset_samples(Temporal,Fraction=='LF')

Spatial<- subset_samples(PSEQ_last_rarefied,!Typeofstudy=='Temporal')
MeuseHW<- subset_samples(Spatial,Typeofstudy=='Spatial_HW') # HW = headwaters
MeuseHW_SP<- subset_samples(MeuseHW,Season=='Spring')
MeuseHW_SU<- subset_samples(MeuseHW,Season=='Summer')
MeuseHW_SP_SF<- subset_samples(MeuseHW_SP,Fraction=='SF')
MeuseHW_SP_LF<- subset_samples(MeuseHW_SP,Fraction=='LF')
MeuseHW_SU_SF<- subset_samples(MeuseHW_SU,Fraction=='SF')
MeuseHW_SU_LF<- subset_samples(MeuseHW_SU,Fraction=='LF')


Spatial_MR<- subset_samples(Spatial,!Typeofstudy=='Spatial_HW') # MR= Main River (so without HW)
Spatial_MR_SP<- subset_samples(Spatial_MR,Season=='Spring')
Spatial_MR_SP_SF<- subset_samples(Spatial_MR_SP,Fraction=='SF')
Spatial_MR_SP_LF<- subset_samples(Spatial_MR_SP,Fraction=='LF')
Spatial_MR_SU<- subset_samples(Spatial_MR,Season=='Summer')
Spatial_MR_SU_SF<- subset_samples(Spatial_MR_SU,Fraction=='SF')
Spatial_MR_SU_LF<- subset_samples(Spatial_MR_SU,Fraction=='LF')


# Alpha diversity measure#### 


  ## Figure 2A - Boxplot alpha diversity headwaters (HW) ####
a_div <- estimate_richness(MeuseHW, measures=c("Shannon"))
# 2)make a dataframe for environmental data
env_data <- as.data.frame(sample_data(MeuseHW))
#combine both datasets( diversity & environmental data)
div_env <- cbind(a_div, env_data)
#Create a new column in metadata combining both fraction size and season
div_envnew<-div_env %>% mutate(Fractionandseason = case_when(Season=="Spring" & Fraction=="SF" ~"SF_spring", Season=="Spring" & Fraction=="LF" ~"LF_spring", Season=="Summer" & Fraction=="SF" ~"SF_summer",Season=="Summer" & Fraction=="LF" ~"LF_summer"))
div_envnew

# Impose the order of the boxplots
env_data$Fractionandseason <- factor(env_data$Fractionandseason, levels = c("Spring_SF","Spring_BF", "Summer_SF","Summer_BF"))


# Add sample size on boxplots
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}


p<-qplot(Fractionandseason, Shannon, data=div_envnew, geom = "boxplot",
         binwidth=0.5, fill=Fraction, main=" Alpha diversity (Shannon index) at the ASV level in headwaters")+ ylab("Alpha diversity (Shannon index)") + geom_boxplot()+ stat_summary(fun.data = give.n, geom = "text", size = 4, position = position_dodge( width = 0.75))+ theme_classic()
p

# Add wilcoxon test comparison two by two
my_comparisons<- list(c("Spring_SF","Spring_LF"),c("Spring_SF", "Summer_SF"),c("Spring_SF", "Summer_LF"),c("Spring_LF", "Summer_SF"), c("Spring_LF", "Summer_LF"), c("Summer_SF", "Summer_LF"))
q<-p + stat_compare_means(label = "p.signif",method = "wilcox.test", comparisons = my_comparisons)
q


  ## Figure 2B - Boxplot CM relative abundance - HW ####

# First, define the Core Microbiome (CM) for each group (Spring SF, Spring LF, Summer SF, Summer LF)
library(microbiome)

# CM Headwater spring SF
s.ra.3<-transform_sample_counts(MeuseHW_SP_SF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_HW_SP_SF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_HW_SP_SF
core_HW_SP_SF@tax_table
df<-data.frame(otu_table (core_TB_SP_F))
write_xlsx(df,"TBSPSF_coremicrobiomeotutable_ASVs")

# CM Headwater spring LF
s.ra.3<-transform_sample_counts(MeuseHW_SP_LF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_HW_SP_LF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_HW_SP_LF
core_HW_SP_LF@tax_table

# CM Headwater summer SF
s.ra.3<-transform_sample_counts(MeuseHW_SU_SF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_HW_SU_SF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_HW_SU_SF
core_HW_SU_SF@tax_table


# CM Headwater summer LF
s.ra.3<-transform_sample_counts(MeuseHW_SU_LF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_HW_SU_LF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_HW_SU_LF
core_HW_SU_LF@tax_table

# Export the otu_table of the 4 HW CM in order to calculate the sum of the relative abundance of all CM ASVs for each sample
df<-data.frame(otu_table (core_HW_SP_SF))
write_xlsx(df,"HW_SP_SF_coremicrobiomeotutable_ASVs")

df<-data.frame(otu_table (core_HW_SP_LF))
write_xlsx(df,"HW_SP_LF_coremicrobiomeotutable_ASVs")

df<-data.frame(otu_table (core_HW_SU_SF))
write_xlsx(df,"HW_SU_SF_coremicrobiomeotutable_ASVs")

df<-data.frame(otu_table (core_HW_SU_LF))
write_xlsx(df,"HW_SU_LF_coremicrobiomeotutable_ASVs")

# Then create and excel file with three columns: the first with the samples names, the second with the sum of the relative abundance of the CM ASVs for each sample (called "Percentage_CM"), and a third being the "Fractionandseason" column already used previously
# This excel file is called "BoxplotHW"

# Then make boxplot

env_data <- sample_data(read.delim2("BoxplotHW.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)

# Impose the order of the boxplots
env_data$Fractionandseason <- factor(env_data$Fractionandseason, levels = c("Spring_SF","Spring_BF", "Summer_SF","Summer_BF"))


##Add p-value to boxplot####
#install.packages("ggpubr")
library(ggpubr)
library (dplyr)

# Add sample size on boxplot
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

p<-qplot(Fractionandseason, Percentage_CM, data=div_envnew, geom = "boxplot", fill=Fractionandseason,
         binwidth=0.5,ylim=c(0,0.8), main="Relative abundance of CM ASVs in HW")+ ylab("Coremicrobiome ASVs percentage") + geom_boxplot()+ stat_summary(fun.data = give.n, geom = "text", size = 4, position = position_dodge( width = 0.75)) +theme_classic()
p

# Add wilcoxon test comparison
my_comparisons<- list(c("Spring_SF","Spring_LF"),c("Spring_SF", "Summer_SF"),c("Spring_SF", "Summer_LF"),c("Spring_LF", "Summer_SF"), c("Spring_LF", "Summer_LF"), c("Summer_SF", "Summer_LF"))
q<-p + stat_compare_means(label = "p.signif",method = "wilcox.test", comparisons = my_comparisons)+scale_fill_manual(values=c("#99CCFF","#0033FF","#FF9999","#CC3300"))
q




library(ggpubr)
library (dplyr)

# Add sample size on boxplot
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

p<-qplot(Fractionandseason, Percent_coremicrobiome, data=env_data, geom = "boxplot", fill=Fractionandseason,
         binwidth=0.5,ylim=c(0,0.8), main="% or core microbiome in headwaters vs Fraction and season")+ ylab("Coremicrobiome ASVs percentage") + geom_boxplot()+ stat_summary(fun.data = give.n, geom = "text", size = 4, position = position_dodge( width = 0.75)) +theme_classic()
p

#  Add p-value
q<-p + stat_compare_means()+ stat_summary(fun.data = give.n, geom = "text", size = 4, position = position_dodge( width = 0.75))
q

# Add wilcoxon test comparison
my_comparisons<- list(c("Spring_SF","Spring_LF"),c("Spring_SF", "Summer_SF"),c("Spring_SF", "Summer_LF"),c("Spring_LF", "Summer_SF"), c("Spring_LF", "Summer_LF"), c("Summer_SF", "Summer_LF"))
q<-p + stat_compare_means(label = "p.signif",method = "wilcox.test", comparisons = my_comparisons)+scale_fill_manual(values=c("#99CCFF","#0033FF","#FF9999","#CC3300"))
q


# Figure 3A####
# Longitudinal profiles of Alpha diversity (Shannon index) at the ASV level - Spatial study MR (Meuse River axis, headwaters excluded) 
a_div <- estimate_richness(Spatial_MR, measures=c("Shannon"))
# 2)make a dataframe for environmental data
env_data <- as.data.frame(sample_data(Spatial_MR))
#combine both datasets( diversity & environmental data)
div_env <- cbind(a_div, env_data)

sp <- ggscatter(div_env, x = "Distance", y = "Shannon",
                 add = "reg.line",               # Add regression line
                 color = "Fractionandseason", palette = c("#99CCFF", "#FF9999","#0033FF", "#CC3300"),
                 shape = "Fraction" , title ='Spatial evolution of the Shannon index on the Meuse river' )   +
  stat_cor(aes(color = Fractionandseason,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.5)+ylim(2.8,6.5) # add correlation equation and set the limits of axes
sp



# Figure 3B######
# Longitudinal profile of the CM relative abundance - Spatial study MR 

# Calculate CM % for the four MR groups


# Spatial MR Spring SF
s.ra.3<-transform_sample_counts(Spatial_MR_SP_SF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_Spatial_MR_SP_SF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_Spatial_MR_SP_SF
core_Spatial_MR_SP_SF@tax_table
df<-data.frame(otu_table (core_Spatial_MR_SP_SF))
write_xlsx(df,"Spatial_MR_SP_SF_coremicrobiomeotutable_ASVs")


# Spatial MR Spring LF
s.ra.3<-transform_sample_counts(Spatial_MR_SP_LF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_Spatial_MR_SP_LF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_Spatial_MR_SP_LF
core_Spatial_MR_SP_LF@tax_table
df<-data.frame(otu_table (core_Spatial_MR_SP_LF))
write_xlsx(df,"Spatial_MR_SP_LF_coremicrobiomeotutable_ASVs")

# Spatial MR Summer SF
s.ra.3<-transform_sample_counts(Spatial_MR_SU_SF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_Spatial_MR_SU_SF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_Spatial_MR_SU_SF
core_Spatial_MR_SU_SF@tax_table
df<-data.frame(otu_table (core_Spatial_MR_SU_SF))
write_xlsx(df,"Spatial_MR_SU_SF_coremicrobiomeotutable_ASVs")


# Spatial MR Summer LF
s.ra.3<-transform_sample_counts(Spatial_MR_SU_LF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_Spatial_MR_SU_LF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_Spatial_MR_SU_LF
core_Spatial_MR_SU_LF@tax_table
df<-data.frame(otu_table (core_Spatial_MR_SU_LF))
write_xlsx(df,"Spatial_MR_SU_LF_coremicrobiomeotutable_ASVs")

# Then create an excel file with 6 columns: Column1= Sample name, C2= distance, C3= relative abundance (in %) of CM "CM_percentage", C4= fraction size, C5= season, C6=Fractionandseason
# We call it CM_Spatial_MR
div_env <- sample_data(read.delim2("CM_Spatial_MR.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) 

library(ggpubr)
sp<- ggscatter(div_env, x = "Distance", y = "CM_percentage",
                 add = "reg.line",               # Add regression line
                 conf.int = TRUE,                # Add confidence interval
                 color = "Fractionandseason", palette = c("#99CCFF", "#FF9999","#0033FF", "#CC3300"),
                 shape = "Fraction" , title ='Longitudinal evolution of the % of core microbiome reads'
)+stat_cor(aes(color = Fractionandseason,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
sp

# Figure 3C#####
# Longitudinal profile of the river discharge - Spatial study MR 

library(ggpubr)
env_data <- as.data.frame(sample_data(Spatial_MR))
sp <- ggscatter(env_data, x = "Distance", y = "Discharge",
                add = "reg.line",               # Add regression line
                conf.int = TRUE,                # Add confidence interval
                color = "Season", palette = c( "#0033FF", "#CC3300"),
                title ='Spatial evolution of the river discharge - MR') 
sp



# Figure 4A####
# Annual profile of Alpha diversity (Shannon index) at the ASV level - Temporal study at Jambes
a_div <- estimate_richness(Temporal, measures=c("Shannon"))
# 2)make a dataframe for environmental data
env_data <- as.data.frame(sample_data(Temporal))
#combine both datasets( diversity & environmental data)
div_env <- cbind(a_div, env_data)

sp <- ggscatter(div_env, x = "Day", y = "Shannon",
                add = "reg.line",               # Add regression line
                color = "Fractionandseason", palette = c("#99CCFF", "#FF9999","#0033FF", "#CC3300"),
                shape = "Fraction" , title ='Temporal evolution of the Shannon index at Jambes - Temporal study' )   +
  stat_cor(aes(color = Fractionandseason,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0.5)+ylim(2.8,6.5) # add correlation equation and set the limits of axes
sp


# Figure 4B####
# Annual profile of the CM relative abundance - Temporal study at Jambes

# Calculate CM % for the 2 fractions


# Temporal SF
s.ra.3<-transform_sample_counts(Temporal_SF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_Temporal_SF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_Temporal_SF
core_Temporal_SF@tax_table
df<-data.frame(otu_table (core_Temporal_SF))
write_xlsx(df,"Temporal_SF_coremicrobiomeotutable_ASVs")

# Temporal LF
s.ra.3<-transform_sample_counts(Temporal_LF, function(x){x / sum(x)}) #transform your asv counts in rel abundances
s.ra.3
core_Temporal_LF<- core(s.ra.3, detection = 0.001, prevalence = 0.9) # here, it was decided to set the detection threshold at 0.1% and the prevalence at 90%
core_Temporal_LF
core_Temporal_LF@tax_table
df<-data.frame(otu_table (core_Temporal_LF))
write_xlsx(df,"Temporal_LF_coremicrobiomeotutable_ASVs")


# Then create an excel file with 6 columns: Column1= Sample name, C2= distance, C3= relative abundance (in %) of CM "CM_percentage", C4= fraction size, C5= season, C6=Fractionandseason
# We call it CM_Spatial_MR
div_env <- sample_data(read.delim2("CM_Temporal.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) 

library(ggpubr)
sp<- ggscatter(div_env, x = "Day", y = "CM_percentage",
               add = "reg.line",               # Add regression line
               conf.int = TRUE,                # Add confidence interval
               color = "Fractionandseason", palette = c("#99CCFF", "#FF9999","#0033FF", "#CC3300"),
               shape = "Fraction" , title ='Temporal evolution of the % of core microbiome reads at Jambes'
)+stat_cor(aes(color = Fractionandseason,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 3)
sp

# Figure 4C####
# Annual profile of the river discharge - Temporal study at Jambes
library(ggpubr)
env_data <- as.data.frame(sample_data(Temporal))
sp <- ggscatter(env_data, x = "Day", y = "Discharge",
                add = "reg.line",               # Add regression line
                conf.int = TRUE,                # Add confidence interval
                color = "Season", palette = c( "#0033FF", "#CC3300"),
                title ='Temporal evolution of the river discharge - Jambes') 
sp




# Figure 5 - Venn Diagram ####
# Number of ASVs of the core microbiomes shared and not shared between the spatial and temporal studies (headwaters excluded). 
# Specific ASVs are those not shared between minimum 2 groups. 

# the otu_tables of the 6 groups that will be compared have already been exported  (Spatial Spring SF- Spatial Spring LF - Spatial Summer SF - Spatial Summer LF - Temporal SF - Temporal LF)
# the spatial groups do not comprise the headwaters (HW) but only the main river (MR)
# they are called Spatial_MR_SP_SF_coremicrobiomeotutable_ASVs, Spatial_MR_SP_LF_coremicrobiomeotutable_ASVs,...

# for each group, a new file is created copying the values of the former otu_tables, and a column "ASV" is added with the number of the ASV, which is visible with the code below
# core_Spatial_MR_SP_SF@tax_table, core_Spatial_MR_SP_LF@tax_table,...
# an "ASVbis" column  is added too next to the "ASV" column with the same values, necessary for performing the following Venn Diagram code

# The same thing is performed for the 6 groups that will be compared
# They are then in 6 separate files

listA<-(read.delim2("Venn_Spatial_MR_SP_SF.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) ## check my excel file example
A<-listA$OTUbis # I had to create this second OTU bis column so that it works...
length(A)

listB<-(read.delim2("Venn_Spatial_MR_SP_LF.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
B<-listB$OTUbis
length(B)

listC<-(read.delim2("Venn_Spatial_MR_SU_SF.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
C<-listC$OTUbis
length(C)

listD<-(read.delim2("Venn_Spatial_MR_SU_LF.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
D<-listD$OTUbis
length(D)

listE<-(read.delim2("Venn_Temporal_SF.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) ## check my excel file example
E<-listE$OTUbis # I had to create this second OTU bis column so that it works...
length(E)

listG<-(read.delim2("Venn_Temporal_LF.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
G<-listG$OTUbis
length(G)


x <- list(A,B,C,D,E,G)

p<-ggVennDiagram(x, label_alpha = 0, lwd = 2,label="count",
                 lty = 'blank',
                 fill = 'black',category.names = c("SF spring","LF Spring","SF summer","LF summer","SF Temporal","LF Temporal"), label_percent_digit = 1) +ggplot2::scale_fill_gradient(low="white",high = "red")+ # I decided to use a red color gradient depending on the amount of unique ASVs of each group 
  ggtitle("Venn Diagram of the core microbiome ASVs - All studies") +theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values=c("black","black","black","black","black","black"))
p




# Figure 6A####
# Relative abundance of core microbiome ASVs in headwaters (regrouped by genera)

# Based the core microbiome of the different HW previously calculated ("core_HW_SP_SF","core_HW_SP_LF","core_HW_SU_SF","core_HW_SP_SF")
# Create a new excel file "CM_HW_Fig6A_otu_table", in which the first column is the ASV number of each CM ASV, the 4 other columns being the relative abundance of each ASV in the 4 different groups (Spring SF, Spring LF, Summer SF, Summer LF)
# The relative abundance of each ASV for each of the 4 groups is calculated based on a geometric mean of their relative abundance among the samples within each of the fours groups
# This file is the otu_table
# Create another file "CM_HW_Fig6A_tax_table", with the first column being the same (number of ASV), and with the taxa information (tax_table), based in the tax_table of the PSEQ objects ("core_HW_SP_SF","core_HW_SP_LF","core_HW_SU_SF","core_HW_SP_SF")

OTU_CM_HW<-(read.delim2("CM_HW_Fig6A_otu_table.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
TAX_CM_HW<-as.matrix(read.delim("CM_HW_Fig6A_tax_table.txt",row.names=1,header=TRUE))

OTU <- otu_table(OTU_CM_HW, taxa_are_rows = TRUE)
TAX <- tax_table(TAX_CM_HW)
# Create a new metadata file "CM_HW_Fig6A_env_data" with four columns, the first being the 4 groups names (HW_SP_SF, HW_SP_LF, HW_SU_SF, HW_SU_LF)
# second column being the fraction size (SF or LF), third column being the group (SpringSF, SpringLF,SummerSF,SummerLF), and the last column being the season
env_data <- sample_data(read.delim2("CM_HW_Fig6A_env_data.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) 

# Impose the order of the representation on the following graphs
env_data$Season <- factor(env_data$Season, levels = c("Spring","Summer"))
env_data$Type <- factor(env_data$Fraction, levels = c("SF","LF"))
env_data$Group <- factor(env_data$Group, levels = c("SpringSF","SpringLF", "SummerSF","SummerLF"))
env_data
CM_HW <- phyloseq(OTU,TAX,env_data)

# The order of the genera is imposed based on their geometric mean among the 4 groups combined
r<-plot_bar(CM_HW,fill="Genus", title='Headwater  ASVs core microbiome pooled by 
            genera ranked by geometric mean-0,1% - 90%')  +ylim(0,0.08)+
  scale_fill_manual(values= c('Limnohabitans'='#01579B','Rhodoferax'='#336633','Flavobacterium'= '#E91E63', 'Yersinia'='#00CCFF','Simplicispira'= '#003300','Cellvibrio'='#000000','Sphaerotilus'='#2196F3',
                              'Pseudomonas'='#FF5722','Methylotenera'='#9900FF','Undibacterium'='#20B2AA','Aeromonas'='#66CCFF','Comamonadaceae_ASV134'='#FFC0CB','Massilia'='#90EE90',"Verrucomicrobiaceae_ASV164"='#330000'))+
  facet_wrap(~Group, scales='free_x',nrow=1) +theme_classic()
r                   

# New desired order to be coherent with top 20
desired_order <- c("Flavobacterium","Rhodoferax","Limnohabitans", "Simplicispira","Pseudomonas","Methylotenera",'Yersinia',"Cellvibrio","Sphaerotilus","Undibacterium","Aeromonas","Comamonadaceae_ASV134","Massilia","Verrucomicrobiaceae_ASV164")

r$data$Genus <- factor(r$data$Genus, levels = desired_order)
print(r)


# Figure 6B####
# Relative abundance of top 20 most abundant genera - Headwaters

# 1) make your phyloseq subsample
MeuseHW

# 2) Export OTU and TAX tables
library(writexl)
df<-data.frame(otu_table(MeuseHW))
write_xlsx(df,"OTU_tb_pooledseasons")
df<-data.frame(tax_table(MeuseHW))
write_xlsx(df,"TAX_tb_pooledseasons")

# Make sum of Spring SF, Spring LF, Summer SF, Summer LF

#2) Rank ASVs  by abundance (make sum of ASVs, add TAX in OTU table, sort them by abundance of ASVs total sum)

#3) For ASVs without name for genus, add in that column the name of the family and the number of the ASV (do this for only top 20, not necessary to do more)
# Example: Comamonodaceae NA --> Comamonodaceae_ASVxx

#3') Make sum of Spring SF, Spring LF, Summer SF, Summer LF, delete samples columns, add those 4 names in metadata file before reimporting this in PSEQ

#4) Reorganize TAX and OTU table in separate files --> only change is the genus name added to the top 20 when no name
# If possible to do this in R, big time gain, but I could'nt figure it out

#5) Reimport tables in R

OTU<-(read.delim2("OTU_HW_last.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
TAX<-as.matrix(read.delim("TAX_HW_last.txt",row.names=1,header=TRUE))
env_data <- sample_data(read.delim2("envdataforbarplotHW.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)
env_data$Season <- factor(env_data$Season, levels = c("Spring", "Summer"))
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
PSEQ <- phyloseq(OTU,TAX,env_data)

PSEQ

# 6) Tax_glom genus --> merge the ASVs in genera
Tax_merged<-tax_glom(PSEQ, "Genus") 


# 7) Take top 20 genus
topN =20 #fine number of most abundant genuses
#get a vector with the names of the OTUs from the most abundant genuses
most_abundant_genus = sort(taxa_sums(Tax_merged), decreasing=TRUE)[1:topN] 
most_abundant_genus

# 8) Subset your Phyloseq Object based on the 20 most abundant genuses
PSEQ_top20= prune_taxa(names(most_abundant_genus),(Tax_merged))
PSEQ_top20
rank_names(PSEQ_top20)

#9) Transform absolute abundance into relative
top20genus <- names(sort(taxa_sums(PSEQ_top20), decreasing=TRUE))[1:20]
ps.top20genus <- transform_sample_counts(PSEQ_top20, function(OTU) OTU/sum(OTU))
ps.top20genus<- prune_taxa(top20genus, ps.top20genus) # maybe this line is useless

# 10) Export the top 20 to see the order of abundance and impose the order on bargraph afterwards
library(writexl)
df<-data.frame(otu_table(ps.top20genus))
write_xlsx(df,"OTU_HW_pstop20genus")
df<-data.frame(tax_table(ps.top20genus))
write_xlsx(df,"TAX_HW_pstop20genus")

# 11) Export the top 20 to divide by amount of total read
library(writexl)
df<-data.frame(otu_table(PSEQ_top20))
write_xlsx(df,"OTU_HW_PSEQ_top20")
df<-data.frame(tax_table(PSEQ_top20))
write_xlsx(df,"TAX_HW_PSEQ_top20")

#11) plot the bargraph
p<-plot_bar(PSEQ_top20, fill="Genus", title='Relative abundance of the top 20 most abundant genus - Upstream rivers - sorted by abundance')+facet_wrap(~Season, scales='free_x',nrow=1)+  theme_classic()
p


#12) Impose colors
q<-p+ scale_fill_manual(values= c('Flavobacterium'= '#E91E63','Rhodoferax'='#336633','Limnohabitans'='#01579B','Simplicispira'= '#003300','Pseudomonas'='#FF5722','Polynucleobacter'='#5C6BC0',
                                  'Acinetobacter'='#00ACC1','Novosphingobium'='#99CC66','Methylotenera'='#9900FF','Comamonadaceae_ASV5'= '#FBC02D','Yersinia'="#66FFFF",'Aurantimicrobium'='#FFFFCC','Cellvibrio'='#000000','Pseudarcobacter'='#CC33CC','TM7a'='#9999FF','Dechloromonas'='#33CCFF','Pseudorhodobacter'='#00FFCC','Fluviicola'='#33FF00',"Pseudarcicella"='#CC33CC','Sphaerotilus'='#2196F3'))

q

#12) Modify the order of bars by checking the top 20 ranking of abundance (otherwise, they appear in alphabetic order)

desired_order <- c("Flavobacterium","Rhodoferax", "Limnohabitans","Simplicispira","Pseudomonas","Polynucleobacter","Acinetobacter","Novosphingobium","Methylotenera","Comamonadaceae_ASV5","Yersinia","Aurantimicrobium","Cellvibrio","Pseudarcobacter","TM7a","Dechloromonas","Pseudorhodobacter","Fluviicola","Pseudarcicella","Sphaerotilus")
q$data$Genus <- factor(q$data$Genus, levels = desired_order)
print(q)



# Figure 6C####
# Same process than for Fig. 6A


OTU_CM_Mainriver<-(read.delim2("CM_Fig6C_otu_table.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
TAX_CM_Mainrriver<-as.matrix(read.delim("CM_Fig6C_tax_table.txt",row.names=1,header=TRUE))

OTU <- otu_table(OTU_CM_Mainriver, taxa_are_rows = TRUE)
TAX <- tax_table(TAX_CM_Mainrriver)
# Create a new metadata file "CM_Fig6C_env_data" with four columns, the first being the 4 groups names (SP_SF, SP_LF, SU_SF, SU_LF)
# second column being the fraction size (SF or LF), third column being the group (SpringSF, SpringLF,SummerSF,SummerLF), and the last column being the season
env_data <- sample_data(read.delim2("CM_Mainriver_Fig6C_env_data.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) 

# Impose the order of the representation on the following graphs
env_data$Season <- factor(env_data$Season, levels = c("Spring","Summer"))
env_data$Type <- factor(env_data$Fraction, levels = c("SF","LF"))
env_data$Group <- factor(env_data$Group, levels = c("SpringSF","SpringLF", "SummerSF","SummerLF"))
env_data
CM_Mainriver <- phyloseq(OTU,TAX,env_data)

Barplot_CM_MR_SP<-plot_bar(CM_Mainriver,fill="Genus", title='ASVs of each core microbiome pooled by 
            genera ranked by average') +facet_wrap(~Typeetude+Type, scales='free_x',nrow=1) +ylim(0,0.8)+
  scale_fill_manual(values= c('Limnohabitans'='#01579B','Flavobacterium'= '#E91E63','Pseudarcicella'= '#FFFFCC','Sporichthyaceae_ASV6'= '#33CCFF','Comamonadaceae_ASV5'= '#FBC02D','Candidatus Planktophila'=  '#33FFCC','hgcI clade'='#4CAF50',
                              'Sediminibacterium'= '#FB8C00','NS11-12 marine group_ASV15'='#660066','Fluviicola'='#33FF00','Sphingorhabdus'='#FBC02D','Algoriphagus'='#CCFFCC','Armatimonas' = '#E91E63', 'Polynucleobacter'='#5C6BC0', 'Rhodoferax'='#336633',
                              'Methylotenera'='#9900FF','Sporichthyaceae_ASV35'='#CCFFFF','Polaromonas'='#003300','Clade III_ASV13'='#CCCCCC','Acidovorax'='#999900','Simplicispira'= '#003300','Tabrizicola'='#FAFAD2','Comamonadaceae_ASV38'='#6633CC',
                              'Rhodobacteraceae_ASV87'='#00FF99','Steroidobacteraceae_ASV127'='#FF0099','Candidatus Limnoluna'='#CCFFFF','Rhodoluna'='#00CC00','Yersinia'='#00CCFF','Candidatus Aquiluna'='#660033',  'NS9 marine group_ASV58'='#CC66FF',
                              'Candidatus Methylopumilus'='#01579B','Verrucomicrobiaceae_ASV164'='#330000','Pseudomonas'='#9900FF','OM60(NOR5) clade'='#ADFF2F',
                              'Thermomonas'='#B0C4DE','Halioglobus'='#00FF00','Aeromonas'='#66CCFF','Comamonadaceae_ASV134'='#FFC0CB'))+
  theme_classic()+theme(legend.position = "none")
Barplot_CM_MR_SP

desired_order <- c('Limnohabitans','Flavobacterium','Pseudarcicella','Sporichthyaceae_ASV6','Comamonadaceae_ASV5','Candidatus Planktophila','hgcI clade',
                   'Sediminibacterium','NS11-12 marine group_ASV15','Fluviicola','Sphingorhabdus','Algoriphagus','Armatimonas', 'Polynucleobacter', 'Rhodoferax',
                   'Methylotenera','Sporichthyaceae_ASV35','Polaromonas','Clade III_ASV13','Acidovorax','Simplicispira','Tabrizicola','Comamonadaceae_ASV38',
                   'Rhodobacteraceae_ASV87','Steroidobacteraceae_ASV127','Candidatus Limnoluna','Rhodoluna','Yersinia','Candidatus Aquiluna',  'NS9 marine group_ASV58',
                   'Candidatus Methylopumilus','Verrucomicrobiaceae_ASV164','Pseudomonas','OM60(NOR5) clade',
                   'Thermomonas','Halioglobus','Aeromonas','Comamonadaceae_ASV134')
Barplot_CM_MR_SP$data$Genus <- factor(Barplot_CM_MR_SP$data$Genus, levels = desired_order)
print(Barplot_CM_MR_SP)


Spatial_Summer

# Impose colors
Barplot_CM_MR_SU<-plot_bar(Spatial_Summer,fill="Genus", title='ASVs of each core microbiome pooled by 
            genera ranked by average-0,1% - 90%') +facet_wrap(~Typeetude+Type, scales='free_x',nrow=1) +ylim(0,0.8)+
  scale_fill_manual(values= c('Limnohabitans'='#01579B','Flavobacterium'= '#E91E63','Pseudarcicella'= '#FFFFCC','Sporichthyaceae_ASV6'= '#33CCFF','Comamonadaceae_ASV5'= '#FBC02D','Candidatus Planktophila'=  '#33FFCC','hgcI clade'='#4CAF50',
                              'Sediminibacterium'= '#FB8C00','NS11-12 marine group_ASV15'='#660066','Fluviicola'='#33FF00','Sphingorhabdus'='#FBC02D','Algoriphagus'='#CCFFCC','Armatimonas' = '#E91E63', 'Polynucleobacter'='#5C6BC0', 'Rhodoferax'='#336633',
                              'Methylotenera'='#9900FF','Sporichthyaceae_ASV35'='#CCFFFF','Polaromonas'='#003300','Clade III_ASV13'='#CCCCCC','Acidovorax'='#999900','Simplicispira'= '#003300','Tabrizicola'='#FAFAD2','Comamonadaceae_ASV38'='#6633CC',
                              'Rhodobacteraceae_ASV87'='#00FF99','Steroidobacteraceae_ASV127'='#FF0099','Candidatus Limnoluna'='#CCFFFF','Rhodoluna'='#00CC00','Yersinia'='#00CCFF','Candidatus Aquiluna'='#660033',  'NS9 marine group_ASV58'='#CC66FF',
                              'Candidatus Methylopumilus'='#01579B','Verrucomicrobiaceae_ASV164'='#330000','Pseudomonas'='#9900FF','OM60(NOR5) clade'='#ADFF2F',
                              'Thermomonas'='#B0C4DE','Halioglobus'='#00FF00','Aeromonas'='#66CCFF','Comamonadaceae_ASV134'='#FFC0CB'))+
  theme_classic()+
  #+ #font white with no grid
  theme(legend.position = "none")
# remove legend
Barplot_CM_MR_SU

desired_order <- c('Limnohabitans','Flavobacterium','Pseudarcicella','Sporichthyaceae_ASV6','Comamonadaceae_ASV5','Candidatus Planktophila','hgcI clade',
                   'Sediminibacterium','NS11-12 marine group_ASV15','Fluviicola','Sphingorhabdus','Algoriphagus','Armatimonas', 'Polynucleobacter', 'Rhodoferax',
                   'Methylotenera','Sporichthyaceae_ASV35','Polaromonas','Clade III_ASV13','Acidovorax','Simplicispira','Tabrizicola','Comamonadaceae_ASV38',
                   'Rhodobacteraceae_ASV87','Steroidobacteraceae_ASV127','Candidatus Limnoluna','Rhodoluna','Yersinia','Candidatus Aquiluna',  'NS9 marine group_ASV58',
                   'Candidatus Methylopumilus','Verrucomicrobiaceae_ASV164','Pseudomonas','OM60(NOR5) clade',
                   'Thermomonas','Halioglobus','Aeromonas','Comamonadaceae_ASV134')
Barplot_CM_MR_SU$data$Genus <- factor(Barplot_CM_MR_SU$data$Genus, levels = desired_order)
print(Barplot_CM_MR_SU)

# Figure 6D####

# Follow the same steps as for Fig. 6B to organize the data
Spatial_MR_SU

# merge the ASVs in genera
Tax_merged<-tax_glom(Spatial_MR_SU, "Genus") 


# Take top 20 genus
topN =20
most_abundant_genus = sort(taxa_sums(Tax_merged), decreasing=TRUE)[1:topN] 
most_abundant_genus

PSEQ_top20= prune_taxa(names(most_abundant_genus),(Tax_merged))
PSEQ_top20
rank_names(PSEQ_top20)


#9) Transform absolute abundance into relative
top20genus <- names(sort(taxa_sums(PSEQ_top20), decreasing=TRUE))[1:20]
ps.top20genus <- transform_sample_counts(PSEQ_top20, function(OTU) OTU/sum(OTU))
ps.top20genus<- prune_taxa(top20genus, ps.top20genus) # maybe this line is useless



# Export the top 20 to see the order of abundance and impose the order on bargraph afterwards

# plot the bargraph
p<-plot_bar(ps.top20genus, fill="Genus", title='Relative abundance of the top 20 most abundant genera - Meuse main river - sorted by abundance')+facet_wrap(~Season, scales='free_x',nrow=1)
p
# Impose colors
q<-p+ scale_fill_manual(values= c('Flavobacterium'= '#E91E63','Limnohabitans'='#01579B', "hgcI clade"='#4CAF50', 'Cyanobium PCC-6307'= '#FAFAFA', 'Candidatus Planktophila'=  '#33FFCC', 'Sporichthyaceae_ASV6'= '#33CCFF' , 'Pseudarcicella'= '#FFFFCC','Fluviicola'='#33FF00','Comamonadaceae_ASV5'= '#FBC02D','Sediminibacterium'= '#FB8C00','Microcystis PCC-7914'='#795548','Polynucleobacter'='#5C6BC0','Rhodoferax'='#336633', 'Clade III_ASV13'='#CCCCCC','NS11-12 marine group_ASV15'='#660066','Sphingorhabdus'='#E53935','Armatimonas' = '#E91E63','Simplicispira'= '#003300','Methylotenera'='#9900FF','Algoriphagus'='#CCFFCC'))
q
# Impose order based on abundance
desired_order <- c("Flavobacterium","Limnohabitans", "hgcI clade", "Cyanobium PCC-6307", "Candidatus Planktophila", "Sporichthyaceae_ASV6" , "Pseudarcicella","Fluviicola","Comamonadaceae_ASV5","Sediminibacterium","Microcystis PCC-7914","Polynucleobacter","Rhodoferax","Clade III_ASV13","NS11-12 marine group_ASV15","Sphingorhabdus","Armatimonas","Simplicispira","Methylotenera","Algoriphagus")
q$data$Genus <- factor(q$data$Genus, levels = desired_order)
print(q)




# Figure 6E####
# Same process than for Figure 6 A and 6C

OTU_CM_Temporal<-(read.delim2("CM_Fig6E_otu_table.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE))
TAX_CM_Temporal<-as.matrix(read.delim("CM_Fig6E_tax_table.txt",row.names=1,header=TRUE))

OTU <- otu_table(OTU_CM_Temporal, taxa_are_rows = TRUE)
TAX <- tax_table(TAX_CM_Temporal)
# Create a new metadata file "CM_Fig6E_env_data" with four columns, the first being the 4 groups names (SP_SF, SP_LF, SU_SF, SU_LF)
# second column being the fraction size (SF or LF), third column being the group (SpringSF, SpringLF,SummerSF,SummerLF), and the last column being the season
env_data <- sample_data(read.delim2("CM_Temporal_Fig6E_env_data.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) 

# Impose the order of the representation on the following graphs
env_data$Season <- factor(env_data$Season, levels = c("Spring","Summer"))
env_data$Type <- factor(env_data$Fraction, levels = c("SF","LF"))
env_data$Group <- factor(env_data$Group, levels = c("SpringSF","SpringLF", "SummerSF","SummerLF"))
env_data
CM_Temporal <- phyloseq(OTU,TAX,env_data)

# Impose colors
q<-plot_bar(CM_Temporal, fill="Genus", title='TJF of AVs core microbiome pooled by 
            genera ranked by average-0,1% - 90%') +ylim(0,0.6)+
  scale_fill_manual(values= c('Limnohabitans'='#01579B','Flavobacterium'= '#E91E63', 'Pseudarcicella'= '#FFFFCC','Comamonadaceae_ASV5'= '#FBC02D', 'Sporichthyaceae_ASV6'= '#33CCFF','Candidatus Planktophila'=  '#33FFCC', 'NS11-12 marine group_ASV15'='#660066','hgcI clade'='#4CAF50',   
                              'Sediminibacterium'= '#FB8C00','Sphingorhabdus'='#FBC02D','Polaromonas'='#003300','Sporichthyaceae_ASV35'='#CCFFFF'))

q                     

desired_order <- c("Limnohabitans","Flavobacterium", "Pseudarcicella","Comamonadaceae_ASV5","Sporichthyaceae_ASV6","Candidatus Planktophila","NS11-12 marine group_ASV15","hgcI clade","Sediminibacterium","Sphingorhabdus","Polaromonas","Sporichthyaceae_ASV35")
q$data$Genus <- factor(q$data$Genus, levels = desired_order)
print(q)

# Figure 7#####
# RDA representing all samples (headwaters, spatial studies on the main river, and temporal study).  
# (A)biotic parameters measured are presented as vectors. 

# Create a PSEQ object without rarefaction
# For that, use the initial otu table obtained after the DADA2 pipeline, called seq.tab.nochim

OTU<-(read.delim2("OTU_tableseqtabnochim.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) # xlxs file are always transormed in txt file before being loaded in R
TAX<-as.matrix(read.delim("TAX_table.txt",row.names=1,header=TRUE)) # correspond à mon tax_table que j'ai exporté après DADA2 pipeline
env_data <- sample_data(read.delim2("metadata.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)

OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
PSEQ <- phyloseq(OTU,TAX,env_data)
PSEQ

# Remove unwanted taxa ####
PSEQ_unclassified_excluded0= subset_taxa(PSEQ, !Kingdom=='NA')
PSEQ_unclassified_excluded0
PSEQ_unclassified_excluded1= subset_taxa(PSEQ_unclassified_excluded0, !Kingdom=='Eukaryota_unclassified')
PSEQ_unclassified_excluded1
PSEQ_unclassified_excluded2= subset_taxa(PSEQ_unclassified_excluded1, !Kingdom=='Archaea')
PSEQ_unclassified_excluded2
PSEQ_unclassified_excluded3= subset_taxa(PSEQ_unclassified_excluded2, !Order=='Chloroplast')
PSEQ_unclassified_excluded3
PSEQ_unclassified_excluded4= subset_taxa(PSEQ_unclassified_excluded3, !Family=='Mitochondria')

PSEQ_last_notrarefied<-PSEQ_unclassified_excluded4

# Perform log ratio transfomation
colsum.nozero<-(colSums(otu_table(PSEQ_last_notrarefied), na.rm=T) !=0)
colsum.nozero
ps<-as.data.frame(t(otu_table(PSEQ_last_notrarefied)))
nozero<-ps[,colsum.nozero] # matrix of asv with only the ASVs that do not have all zeros 
class(nozero)

#transform the 0s in probability vectors
s1<-zCompositions::cmultRepl(nozero, method = 'CZM', delta = 0.5, output = 'p-counts')
view(s1)


#create an alternative ps object with the asv table with probabilities 
a.prob<-phyloseq(otu_table(s1, taxa_are_rows = FALSE),
                 tax_table(tax_table(PSEQ_last_notrarefied)), taxa_names(taxa_names(PSEQ_last_notrarefied)),
                 sample_data(sample_data(PSEQ_last_notrarefied)))


###Centered log-ratio transformation
PSEQ_last_notrarefied_clr <- microbiome::transform(a.prob, "clr")

# Remove the samples not having value for a certain parameter that wxill be used as vector
PSEQ_not_na <-PSEQ_last_notrarefied_clr %>% subset_samples(!is.na(MES)&!is.na(DCO)&!is.na(Prod)&!is.na(Chloa)&!is.na(Prod))
PSEQ_not_na

bray_not_na <- phyloseq::distance(physeq = PSEQ_not_na, method = "euclidean")
cap_ord <- ordinate(
  physeq = PSEQ_not_na, 
  method = "RDA", 
  distance=bray_not_na,formula = ~Temp+DO+COD+TSM+Production+Chloa, scale=TRUE)

cap_plot <- plot_ordination(
  physeq = PSEQ_not_na, 
  ordination = cap_ord, color="Distancefromthemouth",
  axes = c(1,2),title='RDA of both temporal and spatial campaigns (including headwaters' ) + aes(shape = Season)


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")


# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = RDA1*2, # double by 2 the arrow length
                 yend = RDA2*2, 
                 x = 0, # modify the origins of the arrows
                 y = 0, # idem
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 3 * RDA1, 
                 y = 3 * RDA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)
arrowhead = arrow(length = unit(0.01, "npc"))

# Make a new graphic
x<-cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 3,  
    data = arrowdf,
    show.legend = FALSE)    

x

p<-x + 
  scale_color_gradient(low = "dark blue", high = "light blue")+
  theme_classic()
p


# Figure 8#####
# Spearman correlation matrix between (a)biotic parameters (bold) and the 20 most abundant genera (italic) in the Meuse river spatial study (HW excluded, both seasons aggregated)

# Reuse the rarefied PSEQ object
a_div <- estimate_richness(PSEQ_last_rarefied)

# 2)make a dataframe for environmental data
env_data <- as.data.frame(sample_data(PSEQ_last_rarefied))

#combine both datasets( diversity & environmental data)
div_env <- cbind(a_div, env_data)
view(div_env)
# export div_env
div_envdf = as.data.frame(div_env)
library(writexl)
write_xlsx(div_envdf,"divenv_allsamples_Spearman")


# Merge taxa
Tax_merged<-tax_glom(PSEQ_last_rarefied, "Genus") 

#export tax table to get the names of samples
TAXdf = as(tax_table(Tax_merged), "matrix")
#if(taxa_are_rows(Tax_merged)){TAXtaxmerged <- t(TAXtaxmerged)}
# Coerce to data.frame
TAXdf = as.data.frame(TAXdf)
write_xlsx(TAXdf,"TAX_allsamples_Spearman")

OTUdf = as(otu_table(Tax_merged), "matrix")
# Coerce to data.frame
OTUdf = as.data.frame(OTUdf)
write_xlsx(OTUdf,"OTU_allsamples_Spearman")

# after adding top 20 to env_data, reload
env_data <- sample_data(read.delim2("divenv_allsamples_Spearman.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)


mydatacor<-cor(env_data, use = "complete.obs",method = c("spearman"))
mydatacor

# Show table with p values
library("Hmisc")
mydata.rcorr = rcorr(as.matrix(mydatacor))
mydata.rcorr # show correlation and p value

# Visualizing the correlation matrix with different visual options
library(corrplot)

# p-value matrix calculation
p.mat <- cor.mtest(mydatacor)
head(p.mat)


# p-value calculation

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# Only show significant correlation (p-values<0.01)
p<-corrplot(mydatacor, type="upper", 
            p.mat = p.mat, sig.level = 0.01,order="hclust",
            insig = "blank",tl.col="black", tl.srt=90, diag=FALSE)
p


# Figure S1 A####
# Spearman correlation matrices between physico-chemical parameters and the Shannon diversity index of the SF and LF fractions - Spatial study spring (without HW) 
# Based on a metadata file in which for each sample, there is a column for all parameters measured, plus a column for the Shannon index value of SF, and another with a value of Shannon index of LF

env_data <- sample_data(read.delim2("metadata_FigureS1A.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)


mydatacor<-cor(env_data, use = "complete.obs",method = c("spearman"))
mydatacor

# 7) Show table with p values
library("Hmisc")
mydata.rcorr = rcorr(as.matrix(mydatacor))
mydata.rcorr # show correlation and p value

library(corrplot)

p.mat <- cor.mtest(mydatacor)
head(p.mat)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p<-corrplot(mydatacor, type="upper", 
            p.mat = p.mat, sig.level = 0.05,#order="hclust",
            insig = "blank",tl.col="black", tl.cex=1.4,tl.srt=90, diag=FALSE)
p

# Figure S1 B####
# Spearman correlation matrices between physico-chemical parameters and the Shannon diversity index of the SF and LF fractions -  Spatial study summer (without HW) 

# based on a file organized the same way as for Figure S1A
env_data <- sample_data(read.delim2("metadata_FigureS1B.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)


mydatacor<-cor(env_data, use = "complete.obs",method = c("spearman"))
mydatacor

mydata.rcorr = rcorr(as.matrix(mydatacor))
mydata.rcorr # show correlation and p value

# Matrice de p-value de la corrélation
p.mat <- cor.mtest(mydatacor)
head(p.mat)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p<-corrplot(mydatacor, type="upper", 
            p.mat = p.mat, sig.level = 0.05,order="hclust",
            insig = "blank",tl.col="black", tl.cex=1.4,tl.srt=90, diag=FALSE)
p

# Figure S1 C####
# Spearman correlation matrices between physico-chemical parameters and the Shannon diversity index of the SF and LF fractions -  Temporal stud

# based on a file organized the same way as for Figure S1Alibrary(writexl)
env_data <- sample_data(read.delim2("metadata_FigureS1C.txt", header=TRUE,row.names=1,stringsAsFactors = FALSE)) #use this (why? stringAsFactors function make our data look like numbers, not as factors, which means we can do calculations)

mydatacor<-cor(env_data, use = "complete.obs",method = c("spearman"))
mydatacor

mydata.rcorr = rcorr(as.matrix(mydatacor))
mydata.rcorr 

p.mat <- cor.mtest(mydatacor)
head(p.mat)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p<-corrplot(mydatacor, type="upper", 
            p.mat = p.mat, sig.level = 0.05,#order="hclust",
            insig = "blank",tl.col="black", tl.srt=90, diag=FALSE)
p


# Figure S2 A####
# Relative abundance of the top 20 most abundant genera in the main river during spring 

# See script of Figure 6D
# The only difference is that it is made here of the spring group and not the summer one


# Figure S2 B####
# Relative abundance of the top 20 most abundant genera in the temporal study at Jambes

Tax_merged<-tax_glom(Temporal, "Genus") 


# Take top 20 genus
topN =20 #fine number of most abundant OTUs
#get a vector with the names of the OTUs from the most abundant OTUs
most_abundant_genus = sort(taxa_sums(Tax_merged), decreasing=TRUE)[1:topN] 
most_abundant_genus

#Subset your Phyloseq Object based on most abundant OTUs
PSEQ_top20= prune_taxa(names(most_abundant_genus),(Tax_merged))
PSEQ_top20
rank_names(PSEQ_top20)

# For absolute abundance of reads
p<-plot_bar(PSEQ_top20, fill="Genus", title='Relative abundance of the top 20 most abundant genus - temporal study')+facet_wrap(~Type+Season, scales='free_x',nrow=1)+ 
  scale_fill_manual(values= c('Flavobacterium'= '#E91E63','Limnohabitans'='#01579B','Fluviicola'='#33FF00', "Pseudarcicella"='#FFFFCC','hgcI clade'='#4CAF50','Simplicispira'= '#003300', 'Rhodoferax'='#336633','Candidatus Planktophila'=  '#33FFCC','Comamonadaceae_ASV5'='#FBC02D','Polynucleobacter'='#5C6BC0', 'Sporichthyaceae_ASV6'='#33CCFF','Armatimonas' = '#E91E63','Sediminibacterium'= '#FB8C00','Pseudomonas'='#FF5722','NS11-12 marine group_ASV33'='#FF0066','NS11-12 marine group_ASV15'='#9966CC','Acinetobacter'='#00ACC1','Clade III_ASV13'='#CCCCCC','Planomicrobium'='#006600','Methylotenera'='#9900FF'))
p
# Force the order of bars
desired_order <- c("Flavobacterium","Limnohabitans", "Fluviicola", "Pseudarcicella", "hgcI clade", "Simplicispira" , "Rhodoferax","Candidatus Planktophila","Comamonadaceae_ASV5","Polynucleobacter","Sporichthyaceae_ASV6","Armatimonas","Sediminibacterium","Pseudomonas","NS11-12 marine group_ASV33","NS11-12 marine group_ASV15","Acinetobacter","Clade III_ASV13","Planomicrobium","Methylotenera")
p$data$Genus <- factor(p$data$Genus, levels = desired_order)
print(p)
ggsave(file="output/p.pdf",plot=p, device="pdf",width=15, height=8)

# For relative abundance of top 20
top20genus <- names(sort(taxa_sums(PSEQ_top20), decreasing=TRUE))[1:20]
ps.top20genus <- transform_sample_counts(PSEQ_top20, function(OTU) OTU/sum(OTU))
ps.top20genus<- prune_taxa(top20genus, ps.top20genus)

p<-plot_bar(ps.top20genus, fill="Genus", title='Relative abundance of the top 20 most abundant genus (including NA) Meuse  temporal (Jambes)- sorted by abundance')+facet_wrap(~Type+Season, scales='free_x',nrow=1)+ 
  scale_fill_manual(values= c('Flavobacterium'= '#E91E63','Limnohabitans'='#01579B','Fluviicola'='#33FF00', "Pseudarcicella"='#FFFFCC','hgcI clade'='#4CAF50','Simplicispira'= '#003300', 'Rhodoferax'='#336633','Candidatus Planktophila'=  '#33FFCC','Comamonadaceae_ASV5'='#FBC02D','Polynucleobacter'='#5C6BC0', 'Sporichthyaceae_ASV6'='#33CCFF','Armatimonas' = '#E91E63','Sediminibacterium'= '#FB8C00','Pseudomonas'='#FF5722','NS11-12 marine group_ASV33'='#FF0066','NS11-12 marine group_ASV10'='#9966CC','Acinetobacter'='#00ACC1','Clade III_ASV13'='#CCCCCC','Planomicrobium'='#006600','Methylotenera'='#9900FF'))

# Impose the order of bars
desired_order <- c("Flavobacterium","Limnohabitans", "Fluviicola", "Pseudarcicella", "hgcI clade", "Simplicispira" , "Rhodoferax","Candidatus Planktophila","Comamonadaceae_ASV5","Polynucleobacter","Sporichthyaceae_ASV6","Armatimonas","Sediminibacterium","Pseudomonas","NS11-12 marine group_ASV33","NS11-12 marine group_ASV10","Acinetobacter","Clade III_ASV13","Planomicrobium","Methylotenera")
p$data$Genus <- factor(p$data$Genus, levels = desired_order)
print(p)
ggsave(file="output/p.pdf",plot=p, device="pdf",width=18, height=8)


# Figure S3####

# Same script as Figure 8 but only on the Main river samples
