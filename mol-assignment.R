#setting libraries
library(readr)
library(tidyverse)
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)

#importing the data
proteinGroups <- read_delim("https://raw.githubusercontent.com/moghbaie/L1_CRC_IP_MS/master/Input_data/Fourth_Run_04202019/txt/proteinGroups.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
names(proteinGroups)

#filtering out the contaminants and reverse sequence, selecting the LFQ intensity columns
proteinGroups_clean <- proteinGroups %>%
  filter(!str_detect(`Majority protein IDs`, "CON") &
           !str_detect(`Majority protein IDs`, "REV")) %>%
  subset(select = c(1, 7, 233:262)) 

#set 0 values to NA
proteinGroups_clean[proteinGroups_clean == 0] <- NA

#log transformation
proteinGroups_clean[,3:32] <- log2(proteinGroups_clean[3:32])

#calculating average per condition
avgLFQ <- proteinGroups_clean[,3:32]
n <- 1:ncol(avgLFQ)
ind <- matrix(c(n, rep(NA, 3 - ncol(avgLFQ)%%3)), byrow=TRUE, ncol=3)
ind <- data.frame(t(na.omit(ind)))
avgLFQ <- as.data.frame(do.call(cbind, lapply(ind, function(i) rowMeans(avgLFQ[,i], na.rm=TRUE))))
colnames(avgLFQ) <- c("c10", "c2", "c1", "c9", "c5", 
                      "c4", "c3", "c8", "c7", "c6")

#calculating log2foldchange
#foldch163T <- avgLFQ[,9:10] - avgLFQ[,8]
#foldch159T <- avgLFQ[,6:7] - avgLFQ[,5]
#foldch.all <- cbind(avgLFQ[,5], foldch159T, avgLFQ[,8], foldch163T)
names(avgLFQ)

foldchresult <- as.data.frame(apply(avgLFQ,1, function(x){mean(x[c(1:4,6:7,9:10)]) - mean(x[c(5,8)])}))

#t-test for all records
ttestresult <- apply(foldch.all,1,function(x)try({t.test(x[c(1,4)],x[c(2,3,5,6)])$p.value}, silent = TRUE))
padjustresult <- as.data.frame(p.adjust(ttestresult, method = "bonferroni"))

#consolidating
LFQintensity_clean <- cbind(proteinGroups_clean[,1:2], foldchresult, padjustresult)
colnames(LFQintensity_clean) <- c("Protein.ID", "Gene.Name", "Log2FoldChange", "Adjusted.PValue")

#plotting
# add a column of NAs
LFQintensity_clean$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
LFQintensity_clean$diffexpressed[LFQintensity_clean$Log2FoldChange > 0.6 & LFQintensity_clean$Adjusted.PValue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
LFQintensity_clean$diffexpressed[LFQintensity_clean$Log2FoldChange < -0.6 & LFQintensity_clean$Adjusted.PValue < 0.05] <- "DOWN"

diffplot <- ggplot(data=LFQintensity_clean, aes(x=Log2FoldChange, y=-log10(Adjusted.PValue), col=diffexpressed, label = Gene.Name)) + 
  geom_point() + theme_minimal() + 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("blue", "black", "red")) +
  theme(legend.position = "none")
diffplot

