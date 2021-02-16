### libraries needed
library(readr)
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)

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

# get only the first gene name
GeneName <- data.frame(do.call('rbind', strsplit(as.character(proteinGroups_clean$`Gene names`),';',fixed=TRUE)))
proteinGroups_clean$`Gene names` <- GeneName[,1]

### log2 transformation of the LFQ intensity
# first set 0 values to NA because we get -inf otherwise
proteinGroups_clean[proteinGroups_clean == 0] <- NA
proteinGroups_clean[,3:32] <- log2(proteinGroups_clean[3:32])

# making combination tables of: 
#   tumorIgG vs tumorORF1
#     and tumorORF1 vs normalORF1
#ovary
LFQ144_IgG_ORF1 <- proteinGroups_clean[,1:14]
LFQ144_ORF1_NORF1 <- proteinGroups_clean[,c(1:2,9:14,15:17,24:26)]

#liver
LFQ159_IgG_ORF1 <- proteinGroups_clean[,c(1:2,18:20,21:23)]
LFQ159_ORF1_NORF1 <- proteinGroups_clean[,c(1:2,21:23,15:17)]

#colon
LFQ163_IgG_ORF1 <- proteinGroups_clean[,c(1:2,27:32)]
LFQ163_ORF1_NORF1 <- proteinGroups_clean[,c(1:2,30:32,24:26)]

### calculating average log2 LFQ per condition
# each condition has 3 replicates
avgLFQ <- proteinGroups_clean[,3:32]
n <- 1:ncol(avgLFQ)
ind <- matrix(c(n, rep(NA, 3 - ncol(avgLFQ)%%3)), byrow=TRUE, ncol=3)
ind <- data.frame(t(na.omit(ind)))
avgLFQ <- as.data.frame(do.call(cbind, lapply(ind, function(i) rowMeans(avgLFQ[,i], na.rm=TRUE))))
colnames(avgLFQ) <- c("c10", "c2", "c1", "c9", "c5", 
                      "c4", "c3", "c8", "c7", "c6")
avgLFQ[sapply(avgLFQ, is.nan)] <- NA

# Joining them to the previous tables
LFQ144_IgG_ORF1 <- cbind(LFQ144_IgG_ORF1, avgLFQ[,c(1,2,3,4)])
LFQ144_ORF1_NORF1 <- cbind(LFQ144_ORF1_NORF1, avgLFQ[,c(3,4,5,8)])

LFQ159_IgG_ORF1 <- cbind(LFQ159_IgG_ORF1, avgLFQ[,c(6,7)])
LFQ159_ORF1_NORF1 <- cbind(LFQ159_ORF1_NORF1, avgLFQ[,c(7,5)])

LFQ163_IgG_ORF1 <- cbind(LFQ163_IgG_ORF1, avgLFQ[,c(9,10)])
LFQ163_ORF1_NORF1 <- cbind(LFQ163_ORF1_NORF1, avgLFQ[,c(10,8)])


## calc.log2foldchange of tumour IgG and tumour ORF1
LFQ144_IgG_ORF1$tumor <- rowMeans(LFQ144_IgG_ORF1[,c(15:16)], na.rm = TRUE)
LFQ144_IgG_ORF1$ctrl <- rowMeans(LFQ144_IgG_ORF1[,c(17:18)], na.rm=TRUE)
LFQ144_IgG_ORF1$ctrl[sapply(LFQ144_IgG_ORF1$ctrl, is.nan)] <- 0
LFQ144_IgG_ORF1$log2foldchange <- LFQ144_IgG_ORF1$tumor - LFQ144_IgG_ORF1$ctrl

LFQ159_IgG_ORF1$c3[sapply(LFQ159_IgG_ORF1$c3, is.nan)] <- 0
LFQ159_IgG_ORF1$log2foldchange <- LFQ159_IgG_ORF1$c4 - LFQ159_IgG_ORF1$c3

LFQ163_IgG_ORF1$c6[sapply(LFQ163_IgG_ORF1$c6, is.nan)] <- 0
LFQ163_IgG_ORF1$log2foldchange <- LFQ163_IgG_ORF1$c7 - LFQ163_IgG_ORF1$c6


## calc.log2foldchange of tumour and normal tissue ORF1
LFQ144_ORF1_NORF1$tumor <- rowMeans(LFQ144_ORF1_NORF1[,c(15:16)], na.rm = TRUE)
LFQ144_ORF1_NORF1$ctrl <- rowMeans(LFQ144_ORF1_NORF1[,c(17:18)], na.rm=TRUE)
LFQ144_ORF1_NORF1$ctrl[sapply(LFQ144_ORF1_NORF1$ctrl, is.nan)] <- 0
LFQ144_ORF1_NORF1$log2foldchange <- LFQ144_ORF1_NORF1$tumor - LFQ144_ORF1_NORF1$ctrl

LFQ159_ORF1_NORF1$c5[sapply(LFQ159_ORF1_NORF1$c5, is.nan)] <- 0
LFQ159_ORF1_NORF1$log2foldchange <- LFQ159_ORF1_NORF1$c3 - LFQ159_ORF1_NORF1$c5

LFQ163_ORF1_NORF1$c8[sapply(LFQ163_ORF1_NORF1$c8, is.nan)] <- 0
LFQ163_ORF1_NORF1$log2foldchange <- LFQ163_ORF1_NORF1$c6 - LFQ163_ORF1_NORF1$c8

### computing t-test for each sample

# as there are not enough triplicates the t-test will fail for some samples
# tumour IgG vs tumour ORF
LFQ144_IgG_ORF1$pvalue <- apply(LFQ144_IgG_ORF1,1,function(x)try({t.test(as.numeric(x[3:8]),as.numeric(x[9:14]))$p.value}, silent = TRUE))
LFQ144_IgG_ORF1$padjust <- p.adjust(LFQ144_IgG_ORF1$pvalue, method = "bonferroni")

LFQ159_IgG_ORF1$pvalue <- apply(LFQ159_IgG_ORF1,1,function(x)try({t.test(as.numeric(x[3:5]),as.numeric(x[6:8]))$p.value}, silent = TRUE))
LFQ159_IgG_ORF1$padjust <- p.adjust(LFQ159_IgG_ORF1$pvalue, method = "bonferroni")

LFQ163_IgG_ORF1$pvalue <- apply(LFQ163_IgG_ORF1,1,function(x)try({t.test(as.numeric(x[3:5]),as.numeric(x[6:8]))$p.value}, silent = TRUE))
LFQ163_IgG_ORF1$padjust <- p.adjust(LFQ163_IgG_ORF1$pvalue, method = "bonferroni")


# tumour vs normal tissue
LFQ144_ORF1_NORF1$pvalue <- apply(LFQ144_ORF1_NORF1,1,function(x)try({t.test(as.numeric(x[3:8]),as.numeric(x[9:14]))$p.value}, silent = TRUE))
LFQ144_ORF1_NORF1$padjust <- p.adjust(LFQ144_ORF1_NORF1$pvalue, method = "bonferroni")

LFQ159_ORF1_NORF1$pvalue <- apply(LFQ159_ORF1_NORF1,1,function(x)try({t.test(as.numeric(x[3:5]),as.numeric(x[6:8]))$p.value}, silent = TRUE))
LFQ159_ORF1_NORF1$padjust <- p.adjust(LFQ159_ORF1_NORF1$pvalue, method = "bonferroni")

LFQ163_ORF1_NORF1$pvalue <- apply(LFQ163_ORF1_NORF1,1,function(x)try({t.test(as.numeric(x[3:5]),as.numeric(x[6:8]))$p.value}, silent = TRUE))
LFQ163_ORF1_NORF1$padjust <- p.adjust(LFQ163_ORF1_NORF1$pvalue, method = "bonferroni")


### volcano plots

# making a list of dataframes
listORF <- list(LFQ144_ORF1_NORF1, 
                LFQ159_ORF1_NORF1, 
                LFQ163_ORF1_NORF1)
names(listORF) <- c("144T_ORF1_against_159/163N_ORF1", 
                    "159T_ORF1_against_159N_ORF1", 
                    "163T_ORF1_against 163N_ORF1")

listIgGORF <- list(LFQ144_IgG_ORF1, 
                   LFQ159_IgG_ORF1, 
                   LFQ163_IgG_ORF1)

names(listIgGORF) <- c("144T_IgG_against_144T_ORF1", 
                    "159T_IgG_against_159T_ORF1", 
                    "163T_IgG_against 163T_ORF1")

pdf("PlotsTumor_Ctrl_ORF.pdf")
#pdf("PlotsTumor_IgG_ORF.pdf")
for (i in seq(listORF)){
  
  df = as.data.frame(listORF[[i]])
  
  # establishing which proteins are differentially expressed
  # start with column of NAs
  df$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, "UP" 
  df$diffexpressed[df$log2foldchange > 0.6 & df$padjust < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, "DOWN"
  df$diffexpressed[df$log2foldchange < -0.6 & df$padjust < 0.05] <- "DOWN"
  
  # to see the name of the differentially expressed proteins we can add
  df$label <- NA
  df$label[df$diffexpressed != "NO"] <- df$`Gene names`[df$diffexpressed != "NO"]
  
  # plot
  diffplot <- ggplot(data=df, aes(x=log2foldchange, y=-log10(padjust), col=diffexpressed, label = label)) + 
      geom_point() + theme_minimal() + geom_text() + geom_text_repel() +
      geom_vline(xintercept=c(-0.6, 0.6), col="red") +
      geom_hline(yintercept=-log10(0.05), col="red") +
      scale_color_manual(values=c("blue", "black", "red")) +
      theme(legend.position = "none") +
      xlim(-8, 12) + ylim(0, 5) +
      ggtitle(names(listORF[i]))
    
  print(diffplot)
}
dev.off()


