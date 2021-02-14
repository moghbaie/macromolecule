### libraries needed
library(readr)
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2)

### importing the data
proteinGroups <- read_delim("https://raw.githubusercontent.com/moghbaie/L1_CRC_IP_MS/master/Input_data/Fourth_Run_04202019/txt/proteinGroups.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(proteinGroups) 

### filtering out the contaminants and reverse sequence, selecting the LFQ intensity columns
# used the protein IDs instead of the columns "Reverse" and "Potential contaminant" with "+" selection
# the latter was not working for me
proteinGroups_clean <- proteinGroups %>%
  filter(!str_detect(`Majority protein IDs`, "CON") &
           !str_detect(`Majority protein IDs`, "REV")) %>%
  subset(select = c(1, 7, 233:262)) 

### log2 transformation of the LFQ intensity
# first set 0 values to NA because we get -inf otherwise
proteinGroups_clean[proteinGroups_clean == 0] <- NA

proteinGroups_clean[,3:32] <- log2(proteinGroups_clean[3:32])

### calculating average log2 LFQ per condition
# each condition has 3 replicates
avgLFQ <- proteinGroups_clean[,3:32]
n <- 1:ncol(avgLFQ)
ind <- matrix(c(n, rep(NA, 3 - ncol(avgLFQ)%%3)), byrow=TRUE, ncol=3)
ind <- data.frame(t(na.omit(ind)))
avgLFQ <- as.data.frame(do.call(cbind, lapply(ind, function(i) rowMeans(avgLFQ[,i], na.rm=TRUE))))
colnames(avgLFQ) <- c("c10", "c2", "c1", "c9", "c5", 
                      "c4", "c3", "c8", "c7", "c6")

## calc.log2foldchange 
avgLFQ[sapply(avgLFQ, is.nan)] <- NA

avgLFQ$tumor <- rowMeans(avgLFQ[,c(1:4,6:7,9:10)], na.rm = TRUE) #all tissues are pulled
avgLFQ$ctrl <- rowMeans(avgLFQ[,c(5,8)], na.rm=TRUE)
avgLFQ$ctrl[sapply(avgLFQ$ctrl, is.nan)] <- 0

avgLFQ$log2foldchange <- avgLFQ$tumor - avgLFQ$ctrl

### computing t-test for each sample
# first setting data by different tissues
foldch163T <- avgLFQ[,9:10] - avgLFQ$ctrl
foldch159T <- avgLFQ[,6:7] - avgLFQ$ctrl
foldch144T <- avgLFQ[,1:4] - avgLFQ$ctrl
foldch.tumours <- cbind(foldch144T, foldch159T, foldch163T, avgLFQ[,5], avgLFQ[,8])

# as there are not enough triplicates the t-test will fail for some samples
foldch.tumours$pvalue <- apply(foldch.tumours,1,function(x)try({t.test(x[c(9,10)],x[c(1:8)])$p.value}, silent = TRUE))
foldch.tumours$padjust <- p.adjust(foldch.tumours$pvalue, method = "bonferroni")

### making the volcano plot 
# first consolidating results into one clean table
LFQintensity_clean <- cbind(proteinGroups_clean[,1:2], avgLFQ$log2foldchange, foldch.tumours$padjust)
colnames(LFQintensity_clean) <- c("Protein.ID", "Gene.Name", "Log2FoldChange", "Adjusted.PValue")

# establishing which proteins are differentially expressed
# start with column of NAs
LFQintensity_clean$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, "UP" 
LFQintensity_clean$diffexpressed[LFQintensity_clean$Log2FoldChange > 0.6 & LFQintensity_clean$Adjusted.PValue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, "DOWN"
LFQintensity_clean$diffexpressed[LFQintensity_clean$Log2FoldChange < -0.6 & LFQintensity_clean$Adjusted.PValue < 0.05] <- "DOWN"

diffplot <- ggplot(data=LFQintensity_clean, aes(x=Log2FoldChange, y=-log10(Adjusted.PValue), col=diffexpressed)) + 
  geom_point() + theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("blue", "black", "red")) +
  theme(legend.position = "none") +
  xlim(-8, 12) + ylim(0, 10)
diffplot

# if you want to see the name of the differentially expressed proteins we can add this:
LFQintensity_clean$label <- NA
LFQintensity_clean$label[LFQintensity_clean$diffexpressed != "NO"] <- LFQintensity_clean$Gene.Name[LFQintensity_clean$diffexpressed != "NO"]

# plot version 2
library(ggrepel)
diffplot_2 <- ggplot(data=LFQintensity_clean, aes(x=Log2FoldChange, y=-log10(Adjusted.PValue), col=diffexpressed, label = label)) + 
  geom_point() + theme_minimal() + geom_text() + geom_text_repel() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("blue", "black", "red")) +
  theme(legend.position = "none") +
  xlim(-8, 12) + ylim(0, 10)
diffplot_2
