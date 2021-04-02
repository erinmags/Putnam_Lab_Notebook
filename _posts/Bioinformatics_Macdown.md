---
layout: post
title: Bioinformatics Project
date: '2021-03-25'
categories: Protocols
tags: [DNAseq2, DNA Methylation]
---

# Bioinformatics Project
Objective: Quantify gene expression of DNMT1, DNMT3a, and TET1 in response to multiple stressors (temperature and pH) over time for two coral species (*Montipora capitata* and *Pocillopora acuta*)

Here is [My Github Repository for this project](https://github.com/erinmags/Bioinformatics-Project) where my RScripts and RMarkdowns can also be referred to!

Here is [a presentation](https://docs.google.com/presentation/d/1fCZTX7qgj2g9UDM6KFihH2sCLY77UQmPr1yyDGtBPYI/edit#slide=id.gbfe0f19a8d_0_26) I made outlining my project

# Uniprot
[Uniprot](https://www.uniprot.org)

Objective: Uniprot is used to obtain a reference genome

* Use the software Uniprot to look for enzymes DNMT1, DNMT3a, and TET1 in the search bar
	* Filter by reviewd
	* Filter by gene name
* Choosing Proteins
	* DNMT1: I chose the mouse data because it was the longest sequence (1,620) out of all the options
	* DNMT3a: I chose the mouse to keep as much consistency as possible, while still considering the length (908)
	* TET1: I chose the human here because there was over 100 base pairs that differed (2,136)
* Add to Basket
* Download as FASTA file uncompressed
* Open in Text Edit or in Terminal using less or nano function to view it

# Blast
Objective: Blast is used to find regions of similarity between biological sequences

* Log into Bluewaves
`ssh -l username`
* Navigate to Putnam Lab home directory
* 
```
$ cd /data/putnamlab
$ ls (lists all items in directory)
$ mkdir username (makes directory for yourself)
$ cd usernmae (changes directory to work in your new one)
```

* Make script to run BLAST in Bluewaves

	`$ nano "script_name".sh`

* Insert Header
	* This will give you information on memory, time to run, email, etc.
	* *M. capitata* header:

			
			#!/bin/bash 
			#SBATCH --job-name="Bioinformatics"
			#SBATCH -t 100:00:00
			#SBATCH --export=NONE
			#SBATCH --mai_l-type=BEGIN,END,FAIL
			#SBATCH --mail-user=erinmags@uri.edu
			#SBATCH -D /data/putnamlab/erinmags/
			
			# echo "START" $(date)
			
			# module load BLAST+/2.8.1-foss-2018b
			
			# makeblastdb -in /data/putnamlab/erinmags/Mcapitata_holotranscriptome_data_v1/Montipora_capitata_v1_unalinged_transcripts_coding.pep.faa -input_type fasta -dbtype prot -title Mcap.protein -parse_seqids -out Mcap.protein.blast.db
			
			# blastp -query Enzyme_Names.fasta -db Mcap.protein.blast.db -outfmt 6 -max_target_seqs 5 -evalue 1e-05 -out M.cap.protein.blast.out
			
			# echo "STOP" $(date)
			
	* *P. acuta* header:

		```
			#!/bin/bash
		#SBATCH --job-name="Bioinformatics"
		#SBATCH -t 100:00:00
		#SBATCH --export=NONE
		#SBATCH --mail-type=BEGIN,END,FAIL
		#SBATCH --mail-user=erinmags@uri.edu
		#SBATCH -D /data/putnamlab/erinmags/
		
		echo "START" $(date)
		
		module load BLAST+/2.8.1-foss-2018b
		
		makeblastdb -in /data/putnamlab/erinmags/Pacuta_holotranscriptome_data_v1/Pocillopora_acuta_v1_gene_models_coding.pep.faa  -input_type fasta -dbtype prot -title Pacuta.protein.new -parse_seqids -out Pacuta.protein.new.blast.db
		
		blastp -query Enzyme_Names.fasta -db Pacuta.protein.new.blast.db -outfmt 6 -max_target_seqs 5 -evalue 1e-05 -out Pacuta.protein.new.blast.out
		
		echo "STOP" $(date)
	```

* Next, I created a database from *M. capitata* and *P. acuta* predicted protein fasta files
* This will search all sequences in the reference against the query file on the computer
	* The query is the file from Uniprot that we downloaded before
* We need to double check that the script we made in Terminal is executable
	* To do this, type `$ chmod u+wrx blast.sh`
	* If the script in question is green, it is executable!
* Input sequences (the files from Uniprot)
	`$ nano Enzyme_Names.fasta`
* Now, run the script
	* `$ sbatch blast.sh`
		* You can use the command `$ squeue` to show you the job on Terminal
		* To view the output, use the command `$ less filename`

# Uploading necessary files to bluewaves from local computer
* Open Terminal
* Here, I uploaded the Coral_Holobiont_Data
	* `$ scp -r username@bluewaves.uri.edu: path to file on local computer (space) path on bluewaves`
		* Example: `scp -r erinmags@bluewaves.uri.edu: /Users/erinmagliano/Desktop/Putnam-lab/Bioinformatics/Coral_Holobiont_Data /data/putnamlab/erinmags`
			* Tip: Instead of typing out the entire path from your local computer, you can drag the file of interest into the Terminal window!
* Uploading Metadata file
	* `$ scp -r username@bluewaves.uri.edu: path to file on local computer (space) path on bluewaves`
		* Example: `scp -r erinmags@bluewaves.uri.edu: /Users/erinmagliano/Desktop/Putnam-lab/Bioinformatics/Master_Fragment /data/putnamlab/erinmags`
* Unpacking the files within the Coral Holobiont folder
	* Login to bluewaves
	* cd into the directory where the file is
	* To untar the file:
		* ` $ tar -xzf Coral_Holobiont_Data/Mcapitata_holotranscriptome_data_v1.tar.gz`
		* ` $ tar -xzf Coral_Holobiont_Data/Pacuta_holotranscriptome_data_v!.tar.gz`

# Testing Scripts in R with First 10,000 Counts
Note: In its entirety, the counts file is too big to manipulate locally, so we had to run it on bluewaves. We want to subset the counts file so that we can interactively write and test our script in RStudio (locally) before running it on the full counts file in bluewaves.

Here is [My script for *M. capitata* testing the script with the first 10,000 counts](https://github.com/erinmags/Bioinformatics-Project/blob/main/M.cap_Bioinformatics_R_Markdown.Rmd)

Here is [My script for *P. acuta* testing the script with the first 10,000 counts](https://github.com/erinmags/Bioinformatics-Project/blob/main/Pacuta_Bioinformatics_Markdown.Rmd)

* Open Terminal and navigate to file with read count file in it
	* In order to subset the first 10,000 counts, use the command `head -10000 Mcap_RSEM_readcounts>First_10000_Counts` and create a new file with only the first 10,000 rows from the Mcap_RSEM_readcounts file 
	* Repeat previous step for P.acuta RSEM read counts
		* We are taking the first 10,000 counts here in order to test the script in RStudio, since the gene count matrix as a whole is too big to run on the local computer
	* Use the command `$ ls` to check if the file is there
* Open RStudio
	* First, we need to set up the library
	
	```
	library("RColorBrewer")
	library("ggplot2")
	library("tidyverse")
	library("ggpubr")
	library("DESeq2")
	library("dplyr")
	library("Rmisc")
	library("lmerTest")
	library("lme4")
	library("car")
	```
	* 	We also need to upload and adjust our data for the gene count matrix
	*  Uploading data for *M. capitata*:

	```
	#Uploading the Master Fragment with all the information regarding our samples
	Master_Fragment <- read.csv("Master_Fragment.csv")
	#Importing data and telling it to set gene IDs as a row name instead of column
	First_10000_Counts <- read.delim2("Coral_Holobiont_Data/Mcapitata_holotranscriptome_data_v1/First_10000_Counts",row.names = "Name")
	```
	* We use the function read_delim2() because it the file is in a tab-deliminaed format
		* Here, we ran into a problem - the gene names in the read_count file were different then the gene names in the BLAST file we generated. Therefore, we needed to run BLAST again to generate the correct gene names. The file above is the file we ended up working.
	* Uploading data for *P. acuta*
	
	```
	#Importing Data
	Master_Fragment2 <- read.csv("Master_Fragment.csv") 
	#Importing data and telling it to set gene IDs as a row name instead of column
	First_10000_Counts2 <- read.delim2("Coral_Holobiont_Data/	Pacuta_holotranscriptome_data_v1/First_10000_Counts_Pacuta",row.names = "Name")
	#Levels for Timepoint
	factor(Master_Fragment2$Timepoint,levels=c("0 hour","1 hour","6 hour","12 hour","Day 1","30 hour","Day 2","1 week","2 week","4 week","6 week","8 week","12 week","16 week"))
	#Replacing Column Names
	New.cols <- gsub("Pacuta_[A-Z][A-Z][A-Z][A-Z]_TP[0-9][0-9]_","X",colnames(First_10000_Counts2))
	New.cols <- gsub("Pacuta_[A-Z][A-Z][A-Z][A-Z]_TP[0-9]_","X",New.cols)
	#Adding them back
	colnames(First_10000_Counts2) <- New.cols
	```
	
	* At this point I removed the part of the column name so it was only Fragment ID
		* To do this, I used the gsub() function in R, which removes patterns or replaces patterns with other patterns
	* Adding row names for *M. capitata* 

	```
#Add Row Names
New.rows <- paste(sep = "","X",Master_Fragment$Plug_ID)
rownames(Master_Fragment) <- New.rows
#Removing NA rows for Master Fragment Timepoint because they correspond to bleaching samples
Master_Fragment <- Master_Fragment[!is.na(Master_Fragment$Timepoint), ] 
#Levels for Timepoint
factor(Master_Fragment$Timepoint,levels=c("0 hour","1 hour","6 hour","12 hour","Day 1","30 hour","Day 2","1 week","2 week","4 week","6 week","8 week","12 week","16 week"))
#Rename Dataset as M.cap and put rows in correct order
Mcap.Fragment <- Master_Fragment[colnames(First_10000_Counts),]
```
	
	* Adding row names for *P. acuta*

```
#Add Row Names
New.rows <- paste(sep = "","X",Master_Fragment2$Plug_ID)
rownames(Master_Fragment2) <- New.rows
	
#Rename Dataset as Pacuta and put rows in correct order
Pacuta.Fragment <- Master_Fragment2[colnames(First_10000_Counts2),]
```
	

#Constructing DESeq2 Dataset
Objective: The purpose of DESeq2 is to take the raw counts from RNAseq and normalize them for visualization. 

* *M. capitata* DESeq2 Dataset Construction:

```
Mcap.Fragment$Treatment <- factor(Mcap.Fragment$Treatment, levels = c("ATAC","ATHC","HTAC","HTHC"))
#Check Levels
str(Mcap.Fragment$Treatment) 
#Saving First_10000_Counts with gene names as a different dataset
First_10000_Counts_gene <- rownames(First_10000_Counts)
#Converting Decimals into Integars
First_10000_Counts <- data.frame(lapply(First_10000_Counts, as.numeric), stringsAsFactors = FALSE)
#Set Design of DESeq2 Data
mdds <- DESeqDataSetFromMatrix(countData = round(First_10000_Counts),colData = Mcap.Fragment,design = ~Treatment)
```

* *P. acuta* DESeq2 Dataset Construction:

```
Pacuta.Fragment$Treatment <- factor(Pacuta.Fragment$Treatment, levels = c("ATAC","ATHC","HTAC","HTHC"))
#Check Levels
str(Pacuta.Fragment$Treatment) 
#Saving First_10000_Counts with gene names as a different dataset
First_10000_Counts_gene_p <- rownames(First_10000_Counts2)
#Converting Decimals into Integars
First_10000_Counts2 <- data.frame(lapply(First_10000_Counts2, as.numeric), stringsAsFactors = FALSE)
#Set Design of DESeq2 Data
pdds <- DESeqDataSetFromMatrix(countData = round(First_10000_Counts2),colData = Pacuta.Fragment,design = ~Treatment)
```

* At this point, I tested the DESeq2 dataset in R Studio. However it is a large dataset that will need to be run in Bluewaves.
* The variance stabilizing transformation has an assumption of 4, so we need to check that all the values are less than that or the data will not get normalized.
* *M. capitata* test:

```
SF.mdds <- estimateSizeFactors(mdds)
print(sizeFactors(SF.mdds)) #Check that values are less than 4
#Variance Stabilizing Transformation
gvst_mdds<- vst(mdds, blind=TRUE)
vst_transformed_counts <- data.frame(assay(gvst_mdds))
head(assay(gvst_mdds), 3)
#Adding column to vst_transformed_counts
rownames(vst_transformed_counts) <- First_10000_Counts_gene
```

* *P. acuta* test:

```
SF.pdds <- estimateSizeFactors(pdds)
print(sizeFactors(SF.pdds)) #Check that values are less than 4
#Variance Stabilizing Transformation
gvst_pdds<- vst(pdds, blind=TRUE)
P_vst_transformed_counts <- data.frame(assay(gvst_pdds))
head(assay(gvst_pdds), 3)
#Adding column to vst_transformed_counts
rownames(P_vst_transformed_counts) <- First_10000_Counts_gene_p
```

* At this point I ran the script in Bluewaves for the entire dataset. Remember, the above steps were testing the dataset for the first 10,000 counts!

Here is my [M. cap script for the entire dataset that was done in Bluewaves](https://github.com/erinmags/Bioinformatics-Project/blob/main/M.cap.Bluewaves.Script.R)

Here is my [P. acuta script for the entire dataset that was done in Bluewaves](https://github.com/erinmags/Bioinformatics-Project/blob/main/P.acuta.Bluewaves.Script.R)

#Subsetting for Genes of Interest
* In this step, I merged the variance stabilizing transformed (VST) counts with the output we got from running BLAST.
* For reference, here is the [BLAST output for M. cap](https://github.com/erinmags/Bioinformatics-Project/blob/main/M.cap.protein.new.blast.out) and the [BLAST output for P. acuta](https://github.com/erinmags/Bioinformatics-Project/blob/main/Pacuta.protein.new.blast.out)
* Subsetting for *M. capitata*

```
#Adding column names to the blast output
m.cap.blast <- read.table(file="M.cap.protein.new.blast.out", col.names=c("Enzyme_Name","Gene_ID","pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
as.factor(m.cap.blast$Enzyme_Name)
#Upload new counts matrix
blast_vst_transformed_counts <- read.csv("blast_vst_transformed_counts.csv")
#Removing column
blast_vst_transformed_counts <- select(blast_vst_transformed_counts, -X)
#Extracting gene names of interest- BLUEWAVES
blast_vst_transformed_counts <- as.data.frame(blast_vst_transformed_counts)
```

* Subsetting for *P. acuta*

```
#Adding column names to the blast output
Pacuta.blast.out <- read.table(file="Pacuta.protein.new.blast.out", col.names=c("Enzyme_Name","Gene_ID","pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
as.factor(Pacuta.blast.out$Enzyme_Name)
#Adding gene ID column to gene count matrix
P_vst_transformed_counts$Gene_ID <- rownames(P_vst_transformed_counts)
#Extracting gene names of interest- BLUEWAVES
P_vst_transformed_counts <- as.data.frame(P_vst_transformed_counts)
P_vst_transformed_counts <- P_vst_transformed_counts[which(rownames(P_vst_transformed_counts) %in% Pacuta.blast.out$Gene_ID),]
#Combining blast results with gene count matrix
P_vst_transformed_counts <- merge(P_vst_transformed_counts, Pacuta.blast.out[,c(1,2)], by="Gene_ID")
#Upload new count matrix
blast_pacuta_vst_transformed_counts <- read.csv("blast_pacuta_vst_transformed_counts.csv")
```

#Obtaining Means and Standard Errors
Objective: These will be used in the expression plots

* At this point, the dataset I am using is the genes of interest from the BLAST output
	* Here is the [dataset for *M. capitata*](https://github.com/erinmags/Bioinformatics-Project/blob/main/blast_vst_transformed_counts.csv)
	* Here is the [dataset for *P. acuta*](https://github.com/erinmags/Bioinformatics-Project/blob/main/blast_pacuta_vst_transformed_counts.csv)

* Means and standard errors for *M. capitata*

```
#Obtaining means and standard errors
#Changing the format of the dataframe so we can obtain means and standard deviation
blast_vst_transformed_counts <- pivot_longer(blast_vst_transformed_counts, 2:133, values_to = "exp" )
#Adding new column for Sample_ID
Mcap.Fragment$Sample_ID <- rownames(Mcap.Fragment)
#Changing "name" column to sample ID
colnames(blast_vst_transformed_counts)[3] <- "Sample_ID"
#Subsetting Mcap.Frag for columns I want to merge
Mcap.Fragment_subset <- subset(Mcap.Fragment, select = c("Tank","Treatment","Temperature","CO2","Timepoint","Sample_ID"))
#Chaning "Day 1" in Master Frag to 30 hour
filter(Mcap.Fragment, Plug_ID == "2153")
Mcap.Fragment["X2153", "Timepoint"] <- "30 hour"
filter(Mcap.Fragment, Plug_ID == "2153")
#Merging
blast_vst_transformed_counts <- merge(blast_vst_transformed_counts, Mcap.Fragment_subset, by="Sample_ID")
#Finding average of counts
mean.counts <- summarySE(blast_vst_transformed_counts, measurevar="exp", groupvars=c("Timepoint", "Treatment", "Temperature", "CO2","Enzyme_Name","Gene_ID"))
#Removing Day 1 Timepoint becasue it was replaced with 30 hours
mean.counts <- filter(mean.counts,Timepoint!="Day 1")
#Reordering Timepoints
mean.counts$Timepoint <- factor(mean.counts$Timepoint, levels=c("0 hour","6 hour","12 hour","30 hour","1 week","2 week","4 week","6 week","8 week", "12 week","16 week"))
```

* Means and standard errors for *P. acuta*

```
#Obtaining means and standard deviation
#Removing column
blast_pacuta_vst_transformed_counts <- select(blast_pacuta_vst_transformed_counts, -X)
#Changing the format of the dataframe so we can obtain means and standard deviation
blast_pacuta_vst_transformed_counts <- pivot_longer(blast_pacuta_vst_transformed_counts, 2:120, values_to = "exp" )
#Adding new column for Sample_ID
Pacuta.Fragment$Sample_ID <- rownames(Pacuta.Fragment)
#Changing "name" column to sample ID
colnames(blast_pacuta_vst_transformed_counts)[3] <- "Sample_ID"
#Subsetting Mcap.Frag for columns I want to merge
Pacuta.Fragment_subset <- subset(Pacuta.Fragment, select = c("Tank","Treatment","Temperature","CO2","Timepoint","Sample_ID"))
#Merging
blast_pacuta_vst_transformed_counts <- merge(blast_pacuta_vst_transformed_counts, Pacuta.Fragment_subset, by="Sample_ID")
#Finding average of counts
mean.counts.p <- summarySE(blast_pacuta_vst_transformed_counts, measurevar="exp", groupvars=c("Timepoint", "Treatment", "Temperature", "CO2","Enzyme_Name","Gene_ID"))
#Reordering Timepoints
mean.counts.p$Timepoint <- factor(mean.counts.p$Timepoint, levels=c("0 hour","6 hour","12 hour","30 hour","1 week","2 week","4 week","6 week","8 week", "12 week","16 week"))
```

#Subsetting for the Highest Hits
* Initially after following all the steps above, I made my expression plots for both *M. capitata* and *P. acuta* for DNMT1, DNMT3a, and TET1 for each treatment over time. However, these plots were extremely hard to read since they had a ton of variation. Ultimately, I decided to take just the highest hits for each gene and plot those instead in order to lessen the amount of variation. 
* When looking for the genes with the highest hits, one should take that with the lowest e-value and the highest bitscore. 

* Subsetting for highest *M. cap* hits

```
#DNMT1 M.cap gene with highest hit
Keep_m.cap_1 <- subset(mean.counts,Gene_ID == "Montipora_capitata_TRINITY_DN27146_c0_g1_i6.p1")
#DNMT3a M.cap gene with highest hit
Keep_m.cap_2 <- subset(mean.counts,Gene_ID == "Montipora_capitata_TRINITY_DN265265_c0_g1_i4.p1")
#TET M.cap gene with highest hit
Keep_m.cap_3 <- subset(mean.counts,Gene_ID == "Montipora_capitata_TRINITY_DN423918_c0_g2_i1.p1")
#Merging DNMT1, DNMT3a, and TET M.cap highest hit together
Keep_m.cap <- rbind(Keep_m.cap_1,Keep_m.cap_2,Keep_m.cap_3)
```

* Subsetting for highest *P. acuta* hits

```
#DNMT1 P.acuta gene with highest hit
Keep_p.acuta_1 <- subset(mean.counts.p,Gene_ID == "Pocillopora_acuta_PredGene_TCONS_00045317.p1")
#DNMT3a P.acuta gene with highest hit
Keep_p.acuta_2 <- subset(mean.counts.p,Gene_ID == "Pocillopora_acuta_PredGene_TCONS_00029688.p1")
#TET P.acuta gene with highest hit
Keep_p.acuta_3 <- subset(mean.counts.p,Gene_ID == "Pocillopora_acuta_PredGene_TCONS_00033908.p1")
#Merging DNMT1, DNMT3a, and TET P.acuta highest hit together
Keep_p.acuta <- rbind(Keep_p.acuta_1,Keep_p.acuta_2,Keep_p.acuta_3)
```
#Expression Plots
* Now it is finally time to make our expression plots! Remember, the goal is to see the expression of these genes in response to treatment (temperature and pCO2) over time.
* The four treatments are as follows:
	* ATAC = Ambient Temperature, Ambient CO2
	* ATHC = Ambient Temperature, High CO2
	* HTAC = High Temperature, Ambient CO2
	* HTHC = High Temperature, High CO2

* Plotting *M. capitata*
* To see this plot, click [here!](https://github.com/erinmags/Bioinformatics-Project/blob/main/M.cap.plot.jpeg)

```
#Designating colors for treatments
cols <- c("lightblue", "blue", "pink", "red")
#Plotting
Mcap_plot_hits <- ggplot(Keep_m.cap, aes(x=Timepoint, y = exp, group = Treatment, color = Treatment)) +
  geom_line() + geom_point() +
  facet_grid(~ Enzyme_Name, scales = "free") +
  geom_errorbar(aes(ymin=exp-se, ymax=exp+se), width=.2) +
  scale_color_manual(values = cols) +
  xlab("Timepoint") +
  theme_bw() +
  #ylab(expression(AFDW ~(mg ~cm^-2))) +
  ggtitle("Enzyme Expression M.capitata") +
  theme(legend.position = "left") +
  geom_vline(xintercept = c(9.2), colour="black", linetype="dotted") +
  #theme(plot.title = element_text(face = ‘bold.italic’, size = 12, hjust = 0)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.2, hjust = 1.3));Mcap_plot_hits #Set the text angle

ggsave(file="M.cap.plot.jpeg", Mcap_plot_hits, width = 11, height = 6, units = c("in"))
```

* Plotting *P. acuta*
* * To see this plot, click [here!](https://github.com/erinmags/Bioinformatics-Project/blob/main/P.acuta.plot.jpeg)

```
#Designating colors for treatments
cols <- c("lightblue", "blue", "pink", "red")
#Plotting
Pacuta_hits_plot <- ggplot(Keep_p.acuta, aes(x=Timepoint, y = exp, group = Treatment, color = Treatment)) +
  geom_line() + geom_point() +
  facet_grid(~ Enzyme_Name, scales = "free") +
  geom_errorbar(aes(ymin=exp-se, ymax=exp+se), width=.2) +
  scale_color_manual(values = cols) +
  xlab("Timepoint") +
  theme_bw() +
  #ylab(expression(AFDW ~(mg ~cm^-2))) +
  ggtitle("Enzyme Expression P.acuta") +
  theme(legend.position = "left") +
  geom_vline(xintercept = c(9.2), colour="black", linetype="dotted") +
  #theme(plot.title = element_text(face = ‘bold.italic’, size = 12, hjust = 0)) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.2, hjust = 1.3));Pacuta_hits_plot #Set the text angle

ggsave(file="P.acuta.plot.jpeg", Pacuta_hits_plot, width = 11, height = 6, units = c("in"))
```
	
#Statistical Analysis
* Now, it is time to see if we got any significance within the plots. To do this, I used Linear Mixed Models and Type III ANOVA tests. In order to use these, we have to make sure they meet the assumptions, so I also made Residual Q-Q Plots and tested for Homogeneity of the Variances. 
* Statistical tests for *M. capitata*

```
#New dataset
DNMT1_m <- blast_vst_transformed_counts %>% subset(Gene_ID == "Montipora_capitata_TRINITY_DN27146_c0_g1_i6.p1")
#DNMT1
model_DNMT1_m.cap <- lmer(exp ~ Timepoint*Temperature*CO2 + (1|Tank), na.action=na.omit, data=DNMT1_m)
#Residual test
qqPlot(residuals(model_DNMT1_m.cap)) # qqplot and histogram 
hist(residuals(model_DNMT1_m.cap))
#Levens test
leveneTest(DNMT1_m$exp ~ DNMT1_m$Timepoint*DNMT1_m$Temperature*DNMT1_m$CO2)
#ANOVA
Anova(model_DNMT1_m.cap, ddf="lme4", type='III')
summary(model_DNMT1_m.cap)
#New dataset
DNMT3a_m <- blast_vst_transformed_counts %>% subset(Gene_ID == "Montipora_capitata_TRINITY_DN265265_c0_g1_i4.p1")
#DNMT3a
model_DNMT3a_m.cap <- lmer(exp ~ Timepoint*Temperature*CO2 + (1|Tank), na.action=na.omit, data=DNMT3a_m)
#Residual test
qqPlot(residuals(model_DNMT3a_m.cap)) # qqplot and histogram 
hist(residuals(model_DNMT3a_m.cap))
#Levens test
leveneTest(DNMT3a_m$exp ~ DNMT3a_m$Timepoint*DNMT3a_m$Temperature*DNMT3a_m$CO2)
#ANOVA
Anova(model_DNMT3a_m.cap, ddf="lme4", type='III')
summary(model_DNMT3a_m.cap)
#New dataset
TET1_m <- blast_vst_transformed_counts %>% subset(Gene_ID == "Montipora_capitata_TRINITY_DN423918_c0_g2_i1.p1")
#TET
model_TET_m.cap <- lmer(exp ~ Timepoint*Temperature*CO2 + (1|Tank), na.action=na.omit, data=TET1_m)
#Residual test
qqPlot(residuals(model_TET_m.cap)) # qqplot and histogram 
hist(residuals(model_TET_m.cap))
#Levens test
leveneTest(TET1_m$exp ~ TET1_m$Timepoint*TET1_m$Temperature*TET1_m$CO2)
#ANOVA
Anova(model_TET_m.cap, ddf="lme4", type='III')
summary(model_TET_m.cap)
```

* Statistical tests for *P. acuta*

```
#New dataset
DNMT1_p <- blast_pacuta_vst_transformed_counts %>% subset(Gene_ID == "Pocillopora_acuta_PredGene_TCONS_00045317.p1")
#DNMT1
model_DNMT1_p.acuta <- lmer(exp ~ Timepoint*Temperature*CO2 + (1|Tank), na.action=na.omit, data=DNMT1_p)
#Residual test
qqPlot(residuals(model_DNMT1_p.acuta)) # qqplot and histogram 
hist(residuals(model_DNMT1_p.acuta))
#Levens test
leveneTest(DNMT1_p$exp ~ DNMT1_p$Timepoint*DNMT1_p$Temperature*DNMT1_p$CO2)
#ANOVA
Anova(model_DNMT1_p.acuta, ddf="lme4", type='III')
summary(model_DNMT1_p.acuta)
#New dataset
DNMT3a_p <- blast_pacuta_vst_transformed_counts %>% subset(Gene_ID == "Pocillopora_acuta_PredGene_TCONS_00045317.p1")
#DNMT3a
model_DNMT3a_p.acuta <- lmer(exp ~ Timepoint*Temperature*CO2 + (1|Tank), na.action=na.omit, data=DNMT3a_p)
#Residual test
qqPlot(residuals(model_DNMT3a_p.acuta)) # qqplot and histogram 
hist(residuals(model_DNMT3a_p.acuta))
#Levens test
leveneTest(DNMT3a_p$exp ~ DNMT3a_p$Timepoint*DNMT3a_p$Temperature*DNMT3a_p$CO2)
#ANOVA
Anova(model_DNMT3a_p.acuta, ddf="lme4", type='III')
summary(model_DNMT3a_p.acuta)
#New dataset
TET_p <- blast_pacuta_vst_transformed_counts %>% subset(Gene_ID == "Pocillopora_acuta_PredGene_TCONS_00033908.p1")
#TET1
model_TET_p.acuta <- lmer(exp ~ Timepoint*Temperature*CO2 + (1|Tank), na.action=na.omit, data=TET_p)
#Residual test
qqPlot(residuals(model_TET_p.acuta)) # qqplot and histogram 
hist(residuals(model_TET_p.acuta))
#Levens test
leveneTest(TET_p$exp ~ TET_p$Timepoint*TET_p$Temperature*TET_p$CO2)
#ANOVA
Anova(model_TET_p.acuta, ddf="lme4", type='III')
summary(model_TET_p.acuta)
```
