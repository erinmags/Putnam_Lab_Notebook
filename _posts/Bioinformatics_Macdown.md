---
layout: post
title: Bioinformatics Project
date: '2020-11-30'
categories: Protocols
tags: [DNAseq2]
---

# Bioinformatics Project
Objective: Quantify gene expression of DNMT1, DNMT3a, and TET1 in response to multiple stressors (temperature and pH) over time for two coral species (*Montipora capitata* and *Pocillopora acuta*)

# Step 1: Uniprot
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

# Step 2: Blast
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

* Next, I created a database from *M. capitata* predicted protein fasta file
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
		* ` $ tar -xzf Coral_Holobiont_Data/Mcapitata_holotranscriptome_data_v1.tar.gz
		* ` $ tar -xzf Coral_Holobiont_Data/Pacuta_holotranscriptome_data_v!.tar.gz

# Subset Gene Count Matrix
Note: In its entirety, the counts file is too big to manipulate locally, so we had to run it on bluewaves. We want to subset the counts file so that we can interactively write and test our script in RStudio (locally) before running it on the full counts file in bluewaves.

Here is [My Github Repository for this project](https://github.com/erinmags/Bioinformatics-Project) where my RScripts and RMarkdowns can also be referred to!

* Open Terminal and navigate to file with read count file in it
	* Use the command `$ head Mcap_RSEM_readcounts` to print out the first six rows of the file
	* In order to subset the first 10,000 counts, use the command `head -10000 Mcap_RSEM_readcounts>First_10000_Counts` and create a new file with only the first 10,000 rows from the Mcap_RSEM_readcounts file
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
	```
	* 	We also need to upload and adjust our data for the gene count matrix

	```
	#Uploading the Master Fragment with all the information regarding our samples
	Master_Fragment <- read.csv("Master_Fragment.csv")
	#Importing data and telling it to set gene IDs as a row name instead of column
	First_10000_Counts <- read.delim2("Coral_Holobiont_Data/Mcapitata_holotranscriptome_data_v1/First_10000_Counts",row.names = "Name")
	```
	* We use the function read_delim2() because it the file is in a tab-deliminaed format
		* Here, we ran into a problem - the gene names in the read_count file were different then the gene names in the BLAST file we generated. Therefore, we needed to run BLAST again to generate the correct gene names. 



