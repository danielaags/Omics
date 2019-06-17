## Load bam file using the which argument to ScanBamParam
library(Rsamtools)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(Biostrings)
library(BSgenome)


#Get into the right library
setwd("~/Library/Mobile\ Documents/com~apple~CloudDocs/Documents/Postdoc/Sequencing/miseq/190429/Analysis")

#Read bam files
bamfile <- "des_amplicons.bam"
bf <- BamFile(bamfile)

#Define region of interest, check first .bam header

#des range
which <- IRangesList(pCC1FOS_glyco_td = IRanges(9800, 14000))
#ery range
which <- IRangesList(pCC1FOS_glyco_td = IRanges(23300, 29000))


param <- ScanBamParam(what=scanBamWhat())
#Be less stringent with the max_depth and remove strad ID
p_param <- PileupParam(max_depth=10000, distinguish_strand=FALSE)

#pile-up to get the region of interest
res <- pileup(bf, scanBamParam=param, pileupParam=p_param)

##Special function##
pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.numeric(unique(each$pos))
      tablecounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tablecounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}

#Find the frequencies
res_freq <- pileupFreq(res)
#Select only nt frequency and turn unto tidy data
nt_freq <- select(res_freq, start, A, C, G, T)

#Add sequences
wt_seqs <- readDNAStringSet("ery-des_sequences.txt")
MAGE_seqs <- readDNAStringSet("ery-des_sequences_MAGE.txt")

#Define the region of interest, usually the on covered by the MAGE oligo
wt <- tolower(as.character(wt_seqs[[1]]))
wt <- substring(wt, seq(1, nchar(wt), 1), seq(1,nchar(wt),1))

#Define the MAGE sequence
oligo <- tolower(as.character(MAGE_seqs[[1]]))
oligo <- substring(oligo, seq(1, nchar(oligo), 1), seq(1,nchar(oligo),1))

#eryBVI 23386, 23479
nt_freq <- nt_freq %>% mutate(start = as.numeric(start))
nt_freq_ery <- filter(nt_freq, start > 23386 & start < 23479)


#Function to get the # of reads per site
count_eff <- function(nt_data, wt, oligo) {
#Generate empty data frames to save calculations
base <- 1:nrow(nt_data)
mage <- 1:nrow(nt_data)

#Count the frequency to the wt sequence
for (i in 1:nrow(nt_data)){
  #print(i)
  if (wt[i] == "a"){
    base[i] <- nt_data[i,2]
  } else if (wt[i] == "c"){
    base[i] <- nt_data[i,3]
  } else if (wt[i] == "g"){
    base[i] <- nt_data[i,4]
  } else if (wt[i] == "t"){
    base[i] <- nt_data[i,5]
  }
}
print(base)
#Count the frequency to the MAGE oligo
for (i in 1:nrow(nt_data)){
  #print(j)
  if (oligo[i] == "a"){
    mage[i] <- nt_data[i,2]
  } else if (oligo[i] == "c"){
    mage[i] <- nt_data[i,3]
  } else if (oligo[i] == "g"){
    mage[i] <- nt_data[i,4]
  } else if (oligo[i] == "t"){
    mage[i] <- nt_data[i,5]
  }
}
print(mage)

#Get MAGE efficiency per base
print(mage/base)
eff <- (mage/base)*100
#Integrate in the data frame
nt_data <- cbind(eff, MAGE = oligo[1:nrow(nt_data)], nt_data)
#Form groups to plot
nt_data <- transform(nt_data, 
          group = as.factor(findInterval(eff, seq(0, 110, by=55))))
return(nt_data)
}

nt_eff <- count_eff(nt_freq_ery, wt, oligo)

#Plot
#Create my own breaks
my_breaks <- function(x) { if (max(x) > 50) seq(0, 110, 20) else seq(0, 1, 0.2) }

p <- ggplot(nt_eff, aes(as.integer(start), eff, fill = MAGE)) +
  facet_grid(group~., scales="free") +
  geom_bar(stat = "identity") +
  scale_y_continuous(breaks = my_breaks) 

#Change x-axis to nt sequences and add oligos sequences
p +  scale_x_continuous("sequence", breaks = 1881 : 1975, labels = wt, sec.axis = dup_axis(labels = oligo[1:length(wt)])) +
  labs(title="tylA") 

#p + scale_x_discrete(breaks = seq(1881, 1976, 30), labels = c(1, 30, 60, 90))

#Change the text size
p + theme(axis.text = element_text(angle = , size = 10))
