###################################################################
# EXTRACT CD-HIT RESULTS
###################################################################

library(readr)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

# print file name
print(paste("Working on file: ",args))
filename = substr(args, 1, nchar(args)-1)

# Read only column 1,7,8
cd.s.raw =  read_delim(paste0("./", args[1]), delim = " ",col_types = "ic____ii_", col_names = c("cluster","rep","seqStart", "seqEnd"), progress = F)

# cd.s.raw =  read_delim("/Users/solomon/Desktop/RAD_TE_test/FK16_1201_duffyi.merged.1.120.fa.cd_summary1", delim = " ",col_types = "ic____ii_", col_names = c("cluster","rep","seqStart", "seqEnd"), progress = F)

# remove all rows with just cluster no. This also takes away clusters with just the rep.
cd.s = subset(cd.s.raw, !is.na(seqStart) ) 

# calculate median and depth (+1 because rep is not included here)
cd.s.median = rbind.fill(by(cd.s, cd.s[,"cluster"], function(x) 
  data.frame(cluster = x[1,1], 
             rep = x[1,2],
             seqStart_median = median(x$seqStart), 
             seqEnd_median = median(x$seqEnd), 
             depth = nrow(x)+1)))

# save data
write_delim(cd.s.median, paste0(filename, "2"), col_names = F)



