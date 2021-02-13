#setting libraries
library(readr)
library(tidyverse)
library(stringr)
library(dplyr)
library(reshape2)

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
