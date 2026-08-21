# source activate qiime2-2024.10

cat("Loading libraries")
library(dada2); packageVersion("dada2")
library(data.table); packageVersion("data.table")
library(phyloseq); packageVersion("phyloseq")
library(openssl); packageVersion("openssl")
library(theseus); packageVersion("theseus")
library(tidyverse); packageVersion("tidyverse")
library(Biostrings); packageVersion("Biostrings")
library(ggpubr); packageVersion("ggpubr")
library(ggplot2); packageVersion("ggplot2")
library(biohelper); packageVersion("biohelper")
library(insect); packageVersion("insect")
library(ggthemes); packageVersion("insect")
library(reticulate); packageVersion("reticulate")


# Environment
#cutadapt = "/Users/olivierlaroche/miniconda3/envs/qiime2-2020.11/bin/cutadapt" # path to your cutadapt program. Here I am sourcing it from the qiime miniconda env I have installed
cutadapt = "/home/scientist/miniconda3/envs/qiime2-2024.10/bin/cutadapt" # path to your cutadapt program. Here I am sourcing it from the qiime miniconda env I have installed
blastapp_path = "/mnt/HDD/programs/ncbi-blast-2.16.0+/bin/blastn"
cores=16
path_main <- getwd()
path_fastq <- paste0(path_main,"/fastq_files") # CHANGE ME to the directory containing the fastq files after unzipping.
path_results <- paste0(path_main,"/DADA2_results")
if(!dir.exists(path_results)) dir.create(path_results)

# Target taxa options are 12S, 16S bacteria, rbcL, 18S, COI and ITS (fungi)
#target_taxa = "12S"; forward_primer=c("TTAGATACCCCACTATGC"); reverse_primer=c("TAGAACAGGCTCCTCTAG");minLen = 75 # 12S from Valsecchi et al., 2020 (fish & mammals)
#target_taxa = "16S"; forward_primer=c("GTGYCAGCMGCCGCGGTAA"); reverse_primer=c("CCGYCAATTYMTTTRAGTTT"); tL = c(215,190);minLen = 100 # 16S from Parada et al., 2016, V4-V5
#target_taxa = "16S"; forward_primer=c("CCTACGGGNGGCWGCAG"); reverse_primer=c("GACTACHVGGGTATCTAATCC"); tL = c(226,220);minLen = 100 # 16S (usual fish farms) Klindworth et al. (2013), V3-V4
#target_taxa = "16S"; forward_primer=c("GTGYCAGCMGCCGCGGTAA"); reverse_primer=c("GGGACTACNVGGGTWTCTAAT") # 16S (515F, 806R) Caporaso et al. in 2011 (Illumina primers) V4 and V4-5
#target_taxa = "rbcL"; forward_primer=c("AACAGGTGAAGTTAAAGGTTCATAYTT"); reverse_primer=c("TGTAACCCATAACTAAATCGATCAT"); tL = c(120,100);minLen = 0 # rbcL !!!!! tl may also be 90,90 instead
#target_taxa = "18S"; forward_primer=c("AGGGCAAKYCTGGTGCCAGC"); reverse_primer=c("GRCGGTATCTRATCGYCTT"); tL = c(225,216);minLen = 100 # 18S Zhan
#target_taxa = "18S"; forward_primer=c("AGCAGCCGCGGTAATTCC"); reverse_primer=c("TCCGTCAATTYCTTNAASTTTC"); tL = c(290,270);minLen = 100 # 18S VESPA
#target_taxa = "18S"; forward_primer=c("TTGTACACACCGCCC"); reverse_primer=c("CCTTCYGCAGGTTCACCTAC"); tL = c(225,200);minLen = 50 # 18S V9 (expected length ~ 200bp)
#target_taxa = "18S"; forward_primer=c("TTTGTCTGSTTAATTSCG"); reverse_primer=c("CACAGACCTGTTATTGC"); tL = c(225,200);minLen = 50 # 18S V7 (expected length ~ 200bp)
#target_taxa = "16S_mt"; forward_primer=c("AGACGAGAAGACCCTRTG"); reverse_primer=c("GGATTGCGCTGTTATCCC"); tL = c(NA,NA);minLen = 75 # 16S mitochodrial, Valsecchi et al., 2019
target_taxa = "COI"; forward_primer=c("GGWACWGGWTGAACWGTWTAYCCYCC"); reverse_primer=c("TANACYTCNGGRTGNCCRAARAAYCA"); tL = c(225,216);minLen = 100 # COI
#target_taxa = "ITS"; forward_primer=c("GTGARTCATCGAATCTTTG"); reverse_primer=c("GTCCTCCGCTTATTGATATGC"); tL = c(240,220);minLen = 50 # ITS7-ITS4

c_overlap = "--overlap 15" # Minimum overlap for primer matching with cutadapt

norm = T # Whether the taxonomic assingment needs normalisation or not
sqlFile = "/mnt/HDD/reference_db/ncbi/taxonomy_db/accessionTaxa.sql"
ranks = c("Superkingdom", "Kingdom", "Phylum",  "Class", "Order", "Family", "Genus", "Species")
addExtra = T
clusterASVs = T
phylotree = F
blast_method = "both" # Options are "both", "blastn" and "megablast"


#############################################################################
#############################################################################
#############################################################################
#############################################################################
# 16S
taxo_db_16S = "/mnt/HDD/reference_db/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz"
taxo_db_sp_16S = "/mnt/HDD/reference_db/silva/silva_v138.2_assignSpecies.fa.gz"

# rbcL
# https://github.com/fkeck/diatbarcode
taxo_db_rbcL = "/mnt/HDD/reference_db/diat/diat_barcode_v12_tax_assign_dada2.fa"
taxo_db_sp_rbcL = "/mnt/HDD/reference_db/diat/diat_barcode_v12_sp_assign_dada2.fa"

# 18S
taxo_db_silva_18S = "/mnt/HDD/reference_db/silva/silva_132.18s.99_rep_set.dada2.fa.gz" # path to your reference database
taxo_db_pr2 = "/mnt/HDD/reference_db/pr2/pr2_version_5.0.0_SSU_dada2.fasta.gz" # path to your reference database

# COI 
classifier = "/mnt/HDD/reference_db/insect/COI_classifier.rds"
taxo_db_COI = "/home/john/mkCOInr/Leray_COI/COInr/rdp/trimmed_COInr_rdp_trainseq.fasta"
MetaCOXI_classifier_path = "/mnt/HDD/reference_db/MetaCOXI/MetaCOXI_Seqs_adapted_99_trained.rds"

# 12S 
# No classifier at the moment so solely using blast and megablast with NCBI

# ITS
taxo_db_ITS = "/srv/users/Oli/sh_general_release_dynamic_29.11.2022.fasta"
#############################################################################


# Step 1
# Demultiplexing and primer removal
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fas_Fs_raw <- sort(list.files(path_fastq, pattern="R1_001.fastq.gz", full.names = TRUE))
fas_Rs_raw <- sort(list.files(path_fastq, pattern="R2_001.fastq.gz", full.names = TRUE))

fas_Fs_raw[2]
fas_Rs_raw[2]

# This is our set of primers
FWD <- forward_primer
REV <- reverse_primer
FWD_RC <- dada2:::rc(FWD)
REV_RC <- dada2:::rc(REV)

path_cut <- file.path(path_results, "cutadapt")
if(!dir.exists(path_cut)) dir.create(path_cut)
fas_Fs_cut <- file.path(path_cut, basename(fas_Fs_raw))
fas_Rs_cut <- file.path(path_cut, basename(fas_Rs_raw))

R1_flags <- paste(paste("-g", FWD, collapse = " "), paste("-a", REV_RC, collapse = " "))
R2_flags <- paste(paste("-G", REV, collapse = " "), paste("-A", FWD_RC, collapse = " "))

for(i in seq_along(fas_Fs_raw)) {
  cat("Processing", "-----------", i, "/", length(fas_Fs_raw), "-----------\n")
  system2(cutadapt, args = c(R1_flags, R2_flags,
                             "--discard-untrimmed",
                             "--max-n 0",
                             c_overlap,
                             # Optional strong constraint on expected length
                             #paste0("-m ", 250-nchar(FWD)[1], ":", 250-nchar(REV)[1]),
                             #paste0("-M ", 250-nchar(FWD)[1], ":", 250-nchar(REV)[1]),
                             "-o", fas_Fs_cut[i], "-p", fas_Rs_cut[i],
                             fas_Fs_raw[i], fas_Rs_raw[i]))
}

out_fq1 = ShortRead::countFastq(fas_Fs_raw) %>% dplyr::mutate(width=nucleotides/records)%>%dplyr::select(records,width)%>% rownames_to_column("sample_id") %>% setNames(c("sample_id","input","width"))
out_1 <- out_fq1 %>%dplyr::select(sample_id,input)
out_11 = ShortRead::countFastq(fas_Fs_cut) %>% dplyr::mutate(width=nucleotides/records)%>%dplyr::select(records,width)%>% rownames_to_column("sample_id") %>% setNames(c("sample_id","demultiplexed","width"))
write.table(out_fq1, paste0(path_results,"/inputF.csv"), sep=",", quote = F, row.names = T)
write.table(out_11, paste0(path_results,"/demultiplexedF.csv"), sep=",", quote = F, row.names = T)

#out_fq1_R = ShortRead::countFastq(fas_Rs_raw) %>% dplyr::mutate(width=nucleotides/records)%>%dplyr::select(records,width)%>% rownames_to_column("sample_id") %>% setNames(c("sample_id","input","width"))
out_11_R = ShortRead::countFastq(fas_Rs_cut) %>% dplyr::mutate(width=nucleotides/records)%>%dplyr::select(records,width)%>% rownames_to_column("sample_id") %>% setNames(c("sample_id","demultiplexed","width"))
#write.table(out_fq1_R, paste0(path_results,"/inputR.csv"), sep=",", quote = F, row.names = T)
write.table(out_11_R, paste0(path_results,"/demultiplexedR.csv"), sep=",", quote = F, row.names = T)


#### Inspect read quality profiles ####
fas_Fs_cut_ne=fas_Fs_cut[which(out_11$demultiplexed>100)]
fas_Rs_cut_ne=fas_Rs_cut[which(out_11$demultiplexed>100)]

pF <- plotQualityProfile(sample(fas_Fs_cut_ne, replace = FALSE, size = ifelse(length(fas_Fs_cut_ne) < 100, length(fas_Fs_cut_ne), 100)),aggregate = TRUE) + ggplot2::labs(title = "Forward")
pR <- plotQualityProfile(sample(fas_Rs_cut_ne, replace = FALSE, size = ifelse(length(fas_Rs_cut_ne) < 100, length(fas_Rs_cut_ne), 100)),aggregate = TRUE)+ ggplot2::labs(title = "Reverse")
test = ggarrange(pF,pR, nrow = 2)
ggsave(filename = file.path(path_results, "Read_quality_profile_aggregated.pdf"), plot = test, width = 6, height = 8)


# Step 2: Quality filtering
# Place filtered files in filtered/ subdirectory
cat("\nPerforming quality filtering\n")
path_process <- path_cut # If you skipped primers removal, provide the path to your sequences here
fnFs <- sort(list.files(path_process, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path_process, pattern="R2_001.fastq.gz", full.names = TRUE))
#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S\\d+"), `[`, 1)
filtFs <- file.path(path_results, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_results, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


out_2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(0,0), truncLen=tL,
                       maxN=0, maxEE=c(2,2), truncQ=2, minLen=minLen,
                       rm.phix=TRUE,compress=TRUE, multithread=cores) %>% 
                       as.data.frame() %>% 
                       rownames_to_column("sample_id") %>% 
                       setNames(c("sample_id","demultiplexed","filtered"))
head(out_2)
write.table(out_2, paste0(path_results,"/Qfiltering.csv"), sep=",", quote = F, row.names = T)


# Step 3: Learn error rates
cat("\nLearning error rates\n")
filtFs <- sort(list.files(paste0(path_results, "/filtered"), pattern="_F_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(paste0(path_results, "/filtered"), pattern="_R_filt.fastq.gz", full.names = TRUE))

errF <- learnErrors(filtFs, multithread=cores)
saveRDS(errF, file = paste0(path_results,"/errF.rds"))
errR <- learnErrors(filtRs, multithread=cores)
saveRDS(errR, file = paste0(path_results,"/errR.rds"))

perrF <- plotErrors(errF, nominalQ = TRUE) + ggplot2::labs(title = "Error Forward")
perrR <- plotErrors(errR, nominalQ = TRUE) + ggplot2::labs(title = "Error Reverse")
fig = ggarrange(perrF,perrR, nrow = 2)
ggsave(filename = paste0(path_results,"/Error_rates_learning.pdf"), fig ,width = 6, height = 8)


# Step 4: Dereplication
exists <- file.exists(filtFs)
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- gsub("_F_filt.fastq.gz","", names(derepFs))
names(derepRs) <- gsub("_R_filt.fastq.gz","", names(derepRs))

# Step 5: Sample inference
cat("\nDenoising the data\n")
# Step 5: Sample inference
cat("\nDenoising the data\n")
errF = readRDS(paste0(path_results,"/errF.rds"))
errR = readRDS(paste0(path_results,"/errR.rds"))
dadaFs <- dada(derepFs, err=errF, multithread=cores)
dadaRs <- dada(derepRs, err=errR, multithread=cores)


# Step 6: Merging forward and reverse reads
cat("\nMerging reads\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 0, verbose=TRUE)


# Step 7: Construct feature table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Step 8: Remove chimeras
cat("\nRemoving chimeras\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=cores, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, file = paste0(path_results,"/seqtab.nochim.rds"))


# Inspect distribution of sequence lengths
seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
df = as.data.frame(table(nchar(getSequences(seqtab.nochim))))
df$Var1 = as.numeric(as.character(df$Var1))
# Basic density
p <- ggplot(df, aes(x=Var1)) +
 geom_density()+
  #geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1)+
  ggthemes::theme_clean() +
  ylab("Density") +
  #geom_text(aes(x=mean(length), label=paste0("Mean\n",mean(length)), y=0.25)) +
  scale_x_continuous(name ="Sequence length (bp)",breaks = c(0,50,100,150,200,250,300,350,374,400,450,500,550), labels = c("0","50","100","150","200","250","300","350","374","400","450","500","550"))
cat("\nMean and median read length\n")
mean(df$Var1)
cat("\n")
median(df$Var1)
ggsave(filename =  paste0(path_results,"/Seq_lenght_density.pdf"), plot = p,device = "pdf")



# Step 9: Track reads through pipeline
getN <- function(x) sum(getUniques(x))

out_3 = data.frame(names(dadaFs),sapply(dadaFs, getN)) %>% setNames(c("sample_id","denoised"))
out_4 = data.frame(names(dadaFs),sapply(mergers, getN)) %>% setNames(c("sample_id","merged"))
out_5 = data.frame(rownames(seqtab.nochim),rowSums(seqtab.nochim)) %>% setNames(c("sample_id","nonchim"))

out_1[,'sample_id'] = gsub("(.*)_S\\d+_L001.*","\\1", out_1[,'sample_id'] )
out_11[,'sample_id'] = gsub("(.*)_S\\d+_L001.*","\\1", out_11[,'sample_id'] )
out_2[,'sample_id'] = gsub("(.*)_S\\d+_L001.*","\\1", out_2[,'sample_id'] )

out = dplyr::left_join(out_1,out_2, by="sample_id")
out = dplyr::left_join(out,out_3, by="sample_id")
out = dplyr::left_join(out,out_4, by="sample_id")
out = dplyr::left_join(out,out_5, by="sample_id")
out = out %>% base::replace(is.na(.),0)
out <- out[gtools::mixedorder(out$sample_id, decreasing = T), ]
head(out)

write.table(out, paste0(path_results,"/Read_counts.csv"), sep=",", quote = F, row.names = F)


# Step 10: Saving sequences
cat("\nSaving sequences\n")
uniquesToFasta(getUniques(seqtab.nochim), fout= paste0(path_results,"/uniqueSeqs.fasta"),ids=as.character(as.list(md5(names(getUniques(seqtab.nochim))))))


# Step 11: Assing taxonomy
if(target_taxa == "16S"){
	cat("\nAssigning taxonomy for 16S data with the NB classifier\n")
	seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
	# This uses a Naives-Bayes trained classifier
	taxa <- assignTaxonomy(seqtab.nochim, taxo_db_16S, multithread=cores)
	taxa.plus <- addSpecies(taxa, taxo_db_sp_16S, allowMultiple=TRUE, verbose=TRUE)
	rownames(taxa.plus) <-openssl::md5(rownames(taxa.plus))
	saveRDS(taxa.plus, file = paste0(path_results,"/taxaNB.rds"))
	taxa = taxa.plus
	
	}else if(target_taxa=="ITS"){
	cat("\nAssigning taxonomy for fungi (ITS) data with the NB classifier\n")
	seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
	# This uses a Naives-Bayes trained classifier
	taxa <- assignTaxonomy(seqtab.nochim, taxo_db_ITS, multithread=cores)
	rownames(taxa) <-openssl::md5(rownames(taxa))
	saveRDS(taxa, file = paste0(path_results,"/taxaNB.rds"))
	
	}else if(target_taxa == "rbcL"){
  cat("\nAssigning taxonomy for rbcL data with the NB classifier\n")
  seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
  # This uses a Naives-Bayes trained classifier
  taxa <- assignTaxonomy(seqtab.nochim, taxo_db_rbcL, taxLevels = c("Empire", "Kingdom", "Subkingdom", "Phylum", "Class", "Order", "Family", "Genus"),verbose = TRUE, multithread = cores)
  taxa.plus <- addSpecies(taxa, taxo_db_sp_rbcL, allowMultiple=TRUE, verbose=TRUE)
  rownames(taxa.plus) <-openssl::md5(rownames(taxa.plus))
  saveRDS(taxa.plus, file = paste0(path_results,"/taxaNB.rds"))
  taxa = taxa.plus
  
  }else if(target_taxa == "18S"){
	seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
	cat("\nAssigning taxonomy with NB classifier and Silva 132v\n")
	# This uses a Naives-Bayes trained classifier
	taxa_silva <- assignTaxonomy(seqtab.nochim, taxo_db_silva_18S, multithread=cores)
	rownames(taxa_silva) <-openssl::md5(rownames(taxa_silva))
	taxa_silva <- taxa_silva[,!is.na(colnames(taxa_silva))]
	saveRDS(taxa_silva, file = paste0(path_results,"/taxaSilvaNB.rds"))
	cat("\nAssigning taxonomy with NB classifier and Pr2\n")
	# This uses a Naives-Bayes trained classifier
	taxa_pr2 <- assignTaxonomy(seqtab.nochim, taxo_db_pr2, multithread=cores)
	rownames(taxa_pr2) <-openssl::md5(rownames(taxa_pr2))
	taxa_pr2 <- taxa_pr2[,!is.na(colnames(taxa_pr2))]
	saveRDS(taxa_pr2, file = paste0(path_results,"/taxaPr2NB.rds"))
	
	cat("\nAssigning taxonomy with blast and NCBI/Genbank\n")
	dna_path <- paste0(path_results,"/uniqueSeqs.fasta")
	biohelper::blastn_taxo_assignment(
  		blastapp_path = blastapp_path,
  		method = blast_method,
  		queries = dna_path,
  		db = "nt",
  		output_path = path_results,
  		nthreads = cores,
  		minSim = 100,
  		minCov = 80,
  		update = T,
  		pident = "after",
  		pgenus = 97,
  		pfamily = 95,
  		porder = 87,
  		pclass = 83,
  		pphylum = 81,
  		pkingdom = 79,
  		taxonly="TRUE")
	cat("\nMerging the different assignments\n")
	taxa_silva = readRDS(paste0(path_results,"/taxaSilvaNB.rds")) %>% as.data.frame() %>% tibble::rownames_to_column("feature_id") 
	taxa_pr2 = readRDS(paste0(path_results,"/taxaPr2NB.rds")) %>% as.data.frame() %>% tibble::rownames_to_column("feature_id") 
	blastn = fread(paste0(path_results,"/blastn_taxo_assingment.csv"))
	consensus_taxo = biohelper::taxo_merge(
                      df_list = list(taxa_silva,taxa_pr2,blastn),
                      sqlFile = sqlFile,
                      ranks = c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
                      addExtra = addExtra)
    write.csv(x = consensus_taxo, file = paste0(path_results,"/consensus_taxonomy.csv"),row.names = F)
    taxa = consensus_taxo

}else if(target_taxa %in% c("12S","16S_mt")){	
	cat("\nAssigning taxonomy with blastn and NCBI/Genbank\n")
	dna_path <- paste0(path_results,"/uniqueSeqs.fasta")
	biohelper::blastn_taxo_assignment(
  		blastapp_path = blastapp_path,
  		method = blast_method,
  		queries = dna_path,
  		db = "nt",
  		output_path = path_results,
  		nthreads = cores,
  		minSim = 97,
  		minCov = 80,
  		update = T,
  		pident = "after",
  		pgenus = 95,
  		pfamily = 87,
  		porder = 83,
  		pclass = 81,
  		pphylum = 79,
  		pkingdom = 71,
  		taxonly="TRUE")
	cat("\nMerging the different assignments\n")
	blastn = fread(paste0(path_results,"/blastn_taxo_assingment.csv"))
	#consensus_taxo = biohelper::taxo_merge(
    #                  df_list = list(taxa_silva,taxa_pr2,blastn),
    #                  sqlFile = sqlFile,
    #                  ranks = c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    #                  addExtra = T)
    consensus_taxo=blastn
    colnames(consensus_taxo)[1] = "feature_id"
    write.csv(x = consensus_taxo, file = paste0(path_results,"/consensus_taxonomy.csv"),row.names = F)
    taxa = consensus_taxo
	
}else if(target_taxa == "COI"){
	#cat("\nAssigning taxonomy with NB classifier and COInr_2023_05_03 database (BOLD/NCBI)\n")
	#seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
	#taxa_COInr <- assignTaxonomy(seqtab.nochim, taxo_db_COI, multithread=cores)
	#rownames(taxa_COInr) <-openssl::md5(rownames(taxa_COInr))
	#taxa_COInr <- taxa_silva[,!is.na(colnames(taxa_COInr))]
	#saveRDS(taxa_COInr, file = paste0(path_results,"/taxaCOInrNB.rds"))
	
	cat("\nAssigning taxonomy with DECIPHER MetaCOXI classifier\n")
	MetaCOXI_classifier <- readRDS(MetaCOXI_classifier_path)
	# Read sequences to assign
	seq <- readDNAStringSet(paste0(path_results,"/uniqueSeqs.fasta"))
	# Get the taxonomy from the training set
	ids <- DECIPHER::IdTaxa(seq,MetaCOXI_classifier,type="extended",strand="both", processors=NULL,threshold=0)
	ids
	saveRDS(ids,paste0(path_results,"/ids_MetaCOXI.rds"))
	# Transform to a dataframe
	# Note: ids are provided as a list so we need to transform to dataframe
	n_seq <- length(ids)
	df_rows <- list()
	taxo_levels <- c("root","superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
	max_len <- length(taxo_levels)
	for (i in 1:n_seq) {
  		seq_name <- names(ids[i])
  		taxonomy <- ids[[i]]$taxon
  		confidence <- ids[[i]]$confidence
  		# Extend vectors to the same length
  		taxonomy <- c(taxonomy, rep(NA, max_len - length(taxonomy)))
  		confidence <- c(confidence, rep(NA, max_len - length(confidence)))
  		seq_name <- rep(seq_name, max_len)  # Repeat to match the length
  		df_rows[[i]] <- data.frame(seq_name, taxonomy, confidence, taxo_levels, stringsAsFactors = FALSE)
	}
	dfc <- purrr::reduce(df_rows, bind_rows) %>% 
  		filter(taxo_levels !="root") %>% 
  		pivot_wider(names_from = taxo_levels, values_from = c(taxonomy, confidence))
	saveRDS(dfc,paste0(path_results,"/dfc_MetaCOXI.rds"))
	write.table(dfc,paste0(path_results,"/dfc_MetaCOXI.csv"), quote = F,row.names = F, sep = ",")
	# Define your threshold value
	x <- 0.6  # Change this as needed. !!! NOTE THAT the confidence is the same for all ranks, so if below, assignation for all ranks will be NA
	# Apply the transformation
	dfc <- dfc %>%
  	mutate(across(starts_with("taxonomy_"), 
                ~ ifelse(get(paste0("confidence_", sub("taxonomy_", "", cur_column()))) < x, NA, .),
                .names = "filtered_{.col}"))        
	df_metacoxi <- dfc %>%
  		select(starts_with("taxonomy")) %>%  # Keep only columns starting with "taxonomy"
  		rename_with(~ sub("^taxonomy_", "", .x))  # Remove "taxonomy" from column names
	df_metacoxi=cbind(names(ids),df_metacoxi)
	names(df_metacoxi) = c("asv","superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
	write.table(df_metacoxi,paste0(path_results,"/df_MetaCOXI.csv"), quote = F,row.names = F, sep = ",")

	cat("\nAssigning taxonomy with INSECT COI\n")
	classifier = readRDS(classifier)
	taxa_insect <- classify(colnames(seqtab.nochim), classifier, threshold = 0.8, ping = 0.99, ranks = c("kingdom", "phylum","class","order","family","genus","species"), cores = cores)
	taxa_insect = taxa_insect %>% as.data.frame()
	taxa_insect$representative = md5(colnames(seqtab.nochim))
	write.table(taxa_insect,paste0(path_results,"/taxo_insect.csv"), quote = F,row.names = F, sep = ",")
	saveRDS(taxa_insect, file = paste0(path_results,"/taxid_insect.rds"))
	taxa_insect <- taxa_insect %>% subset(select = c("representative", "kingdom", "phylum","class","order","family","genus","species")) %>% column_to_rownames("representative") %>% as.matrix()
	
	cat("\nAssigning taxonomy with blastn and NCBI/Genbank\n")
	dna_path <- paste0(path_results,"/uniqueSeqs.fasta")
	biohelper::blastn_taxo_assignment(
  		blastapp_path = blastapp_path,
  		method = blast_method,
  		queries = dna_path,
  		db = "nt",
  		output_path = path_results,
  		nthreads = cores,
  		minSim = 97,
  		minCov = 80,
  		update = T,
  		pident = "after",
  		pgenus = 95,
  		pfamily = 87,
  		porder = 83,
  		pclass = 81,
  		pphylum = 79,
  		pkingdom = 71,
  		taxonly="TRUE")
	cat("\nMerging the different assignments\n")
	taxa_insect = fread(paste0(path_results,"/taxo_insect.csv")) %>% dplyr::select(-c(taxID,taxon,rank,score))
	df_metacoxi = fread(paste0(path_results,"/df_MetaCOXI.csv"))
	blastn = fread(paste0(path_results,"/blastn_taxo_assingment.csv"))
	#taxa_COInr = taxa_COInr %>% as.data.frame() %>% tibble::rownames_to_column("feature_id") 
	#taxa_COInr = readRDS(paste0(path_results,"/taxaCOInrNB.rds"))

	consensus_taxo = biohelper::taxo_merge(
                      df_list = list(taxa_insect,blastn,df_metacoxi),
                      sqlFile = sqlFile,
                      ranks = c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
                      addExtra = addExtra)
    write.csv(x = consensus_taxo, file = paste0(path_results,"/consensus_taxonomy.csv"),row.names = F)
    taxa = consensus_taxo
}


# Step 12: Handsoff to Phyloseq
cat("\nCreating a phyloseq object\n")
seqtab.nochim = readRDS(paste0(path_results,"/seqtab.nochim.rds"))
colnames(seqtab.nochim) <-md5(colnames(seqtab.nochim))
taxa_ps = tax_table(taxa %>% as.matrix())
rownames(taxa_ps) = if('feature_id' %in% names(taxa)){taxa[['feature_id']]}else{ rownames(taxa)}
dna <- readDNAStringSet(paste0(path_results,"/uniqueSeqs.fasta")) # Create a DNAStringSet from the ASVs

ps_run <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  taxa_ps,
  refseq(dna))

# Taxo normalisation if needed
if(norm==T){
  ps_run = taxo_normalisation(obj = ps_run, sqlFile = sqlFile, ranks = ranks, addExtra = addExtra)
}

ps_run

otab = pstoveg_otu(ps_run) %>% t() %>% as.data.frame()
otab = cbind(ASVs = rownames(otab), otab)
write.table(otab, paste0(path_results,"/dada_table.txt"),quote=FALSE,sep="\t", row.names = FALSE)

taxa = ps_run@tax_table@.Data %>% as.data.frame()
write.table(taxa, paste0(path_results,"/tax_table.txt"),quote=FALSE,sep="\t", row.names = F)

saveRDS(ps_run, file = paste0(path_results,"/ps_run.rds"))

# Creating a phylogenetic tree
if(phylotree == T){
  cat("\nAligning sequences\n")
  system2("muscle", args = c("-super5",
                           paste0(path_results,"/uniqueSeqs.fasta"),
                           paste0("-output ",path_results,"/aligned_ref_seq.fasta")))
  t = Biostrings::readDNAMultipleAlignment(paste0(path_results,"/aligned_ref_seq.fasta"))
  
  cat("\nComputing distance matrix\n")  
  my_alignment_sequence <- msa::msaConvert(t, type="seqinr::alignment")
  distance_alignment <- seqinr::dist.alignment(my_alignment_sequence)
  
  cat("\nComputing phylogenetic tree using neighbor joining\n")  
  Tree <- ape::bionj(distance_alignment)
  ps_run@phy_tree = phyloseq::phy_tree(Tree)

  cat("\nRooting the tree using the longest branch\n")  
  pick_new_outgroup <- function(tree.unrooted){
    require("magrittr")
    require("data.table")
    require("ape") # ape::Ntip
    # tablify parts of tree that we need.
    treeDT <- 
      cbind(
        data.table(tree.unrooted$edge),
        data.table(length = tree.unrooted$edge.length)
      )[1:Ntip(tree.unrooted)] %>% 
      cbind(data.table(id = tree.unrooted$tip.label))
    # Take the longest terminal branch as outgroup
    new.outgroup <- treeDT[which.max(length)]$id
    return(new.outgroup)
  }
new.outgroup = pick_new_outgroup(ps_run@phy_tree)
ps_run@phy_tree = ape::root(ps_run@phy_tree, outgroup=new.outgroup, resolve.root=TRUE)
saveRDS(ps_run, file = paste0(path_results,"/ps_run.rds"))
}

# Creating OTUs
if(clusterASVs==T){
    cat("\nCreating OTUs\n")
	#dna = ape::read.FASTA(paste0(path_results,"/uniqueSeqs.fasta"))
	dna = ps_run@refseq %>% ape::as.DNAbin()
	
	clusters <- kmer::otu(dna, threshold = 0.97, nstart = 20, k = 5) %>% 
	as.data.frame() %>% rownames_to_column("cluster")
	write_csv(x = clusters, paste0(path_results,"/otu-clusters.csv"))	
	
	clusters$feature_id = gsub("\\*","", clusters$cluster)
	names(clusters)[2] = "otu_id"
	
	clusters = clusters %>%
	dplyr::group_by(otu_id) %>%
	dplyr::mutate(otu_id = cluster[stringi::stri_detect_fixed(cluster, "*")])
	clusters$otu_id = gsub("\\*","", clusters$otu_id)
	cat("\nNumber of ASVs\n")
	print(clusters$feature_id %>% unique() %>% length())
	cat("\nNumber of OTUs\n")
	print(clusters$otu_id %>% unique() %>% length())
	clusters = clusters[match(taxa_names(ps_run),clusters$feature_id),]
	ps_run@tax_table@.Data = ps_run@tax_table@.Data %>% as.data.frame() %>% dplyr::mutate("otu" = clusters$otu_id) %>% as.matrix()
	saveRDS(ps_run, file = paste0(path_results,"/ps_run.rds"))

	# Creating a new otu and table and fasta file	
	# seqtab needs to have rows = ASVs and columns = samples, then merging ASVs in the same cluster, summing abundances within samples
	#seqtab = fread(paste0(path_results,"/dada_table.txt"))
	#names(seqtab)[1] = "feature_id"
	seqtab = ps_run %>% pstoveg_otu() %>% t() %>% as.data.frame() %>% rownames_to_column("feature_id")

	merged_seqtab <- seqtab %>%
  		dplyr::left_join(clusters %>% dplyr::select(-cluster), by = "feature_id") %>%
  		dplyr::group_by(otu_id) %>%
  		dplyr::summarize_at(vars(-feature_id), sum)

	write_csv(x = merged_seqtab, file = paste0(path_results,"/otu-table.csv"))
	dna_t = dna[names(dna) %in% merged_seqtab$otu_id]
	ape::write.FASTA(x = dna_t ,file= paste0(path_results,"/otu-dna-sequences.fasta"))
}


