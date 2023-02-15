library("dada2")
require(stringr)



path <- "C:/Users/svdz/Documents/stage/data/raw_data" 
list.files(path)



###   Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_L001_R1_001_000000000", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001_000000000", full.names = TRUE))



###   remove incorrect samples
fnFs <- fnFs[!grepl("Bean",fnFs)]
fnFs <- fnFs[!grepl("banana",fnFs)]
fnRs <- fnRs[!grepl("Bean",fnRs)]
fnRs <- fnRs[!grepl("banana",fnRs)]



###   Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
###   sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`)
sample.names <- word(basename(fnFs), 1, sep = "_[0-9]{5}")



###   The mean quality score at each position is shown by the green line, 
###   and the quartiles of the quality score distribution by the orange lines.
#plotQualityProfile(fnFs[1:2]) #quality profile forward reads
#plotQualityProfile(fnRs[1:2]) #quality profile backward reads



###   Place filtered files in filtered/ subdirectory 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncLen=c(250,200), trimLeft = c(17,21), 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
#head(out)
#write.csv(out, "out.csv")



###   error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)



###   sample interference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]



###   merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#Inspect the merger data.frame from the first sample
head(mergers[[1]])



###   construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))



###   Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)



###   Track reads trhough pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)



###   Taxonomy
Silva <- "~/tax/Silva/silva_nr99_v138.1_train_set.fa.gz"
RDP <- "~/tax/RDP/rdp_train_set_18.fa.gz"
SilvaS <- "~/tax/Silva/silva_species_assignment_v138.1.fa.gz"
RDPS <- "~/tax/RDP/rdp_species_assignment_18.fa.gz"

taxa <- assignTaxonomy(seqtab.nochim, Silva, multithread=TRUE)
taxa <- addSpecies(taxa, SilvaS)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)










