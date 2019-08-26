#########################################
# RAD_TE_summary
#########################################

args <- commandArgs(trailingOnly = TRUE)

print(paste0("Input file 1: " ,args[1]))
print(paste0("Input file 2: " ,args[2]))
print(paste0("Input file 3: " ,args[3]))

# testing purpose
# args = c("small.fa.cd.summary2", "small.fa.cd.out", "small.fa.cd.masked.annot")


library(readr)
library(plyr)
library(fitdistrplus)

#~~~~~~~~~~~~~~~~~~~~~~~~
print("Read cd-hit's .summary2")
#~~~~~~~~~~~~~~~~~~~~~~~~
# read cd.summary2
cdhit = read_delim(paste0("./",args[1]), delim = " ",col_types = "ncnnn" ,col_names = c("cluster", "rep","seqStart", "seqEnd", "depth"), progress = F)

# change TE.tag to 18 characters
cdhit$rep18 = strtrim(cdhit$rep,18)


#~~~~~~~~~~~~~~~~~~~~~~~~
print("Read RepeatMasker's .cd.out")
#~~~~~~~~~~~~~~~~~~~~~~~~
RMout = read_table2(paste0("./",args[2]), skip=3, col_types = "____cnn__cc____c", col_names = c("TE.tag", "repStart", "repEnd", "repeat", "class_family", "remove"), progress = F)

# Remove all * in remove 
RMout = subset(RMout, is.na(remove))[-6]

# Calculate TE length
RMout$bp = RMout$repEnd - RMout$repStart + 1

# Remove duplicate: some TE.tag have two hits. Use the one with longer coverage
RMout = as.data.frame(RMout)
temp = RMout[order(RMout[,"TE.tag"],-RMout[,"bp"]),]
RMout = temp[!duplicated(temp$TE.tag),]
rm(temp)

# change TE.tag to 18 characters
RMout$TE.tag18 = strtrim(RMout$TE.tag,18)

#~~~~~~~~~~~~~~~~~~~~~~~~
print("Read RepeatProteinMask's .cd.masked.annot")
#~~~~~~~~~~~~~~~~~~~~~~~~
# 
Protein = read_table(paste0("./",args[3]), skip=3, col_types = "?__cnn_cc__", col_names = c("pValue","TE.tag", "repStart", "repEnd", "repeat", "class_family"), progress = F)

# Subet p value to < 0.00001
Protein = subset(Protein, pValue < 0.00001)[-1]

# Calcualte TE length
Protein$bp = Protein$repEnd - Protein$repStart + 1 

# Remove duplicates: some TE.tag have two hits. Use the one with longer coverage
Protein = as.data.frame(Protein)
temp = Protein[order(Protein[,"TE.tag"],-Protein[,"bp"]),]
Protein = temp[!duplicated(temp$TE.tag),]
rm(temp)

# ** TE.tag is trancated: it only has 18 characters. Should be OK--just pay attation when matching
Protein$TE.tag18 = Protein$TE.tag


#~~~~~~~~~~~~~~~~~~~~~~~~
print("Combind data")
#~~~~~~~~~~~~~~~~~~~~~~~~
cdhitTE = merge(cdhit, rbind(RMout, Protein), by.x="rep18", by.y="TE.tag18", all = T)

# The ones that are in TE but not in cdhit are ones with just 1 copy
# Need to keep them for accounting for total library depth!
cdhitTE$depth[is.na(cdhitTE$depth)] = 1

# Make a master tag column
cdhitTE$tag = apply(data.frame(cdhitTE$rep18, cdhitTE$rep, cdhitTE$TE.tag), 1, function(x) x[which.max(nchar(x))])
cdhitTE = cdhitTE[-c(1,3,7)]

# Split class and family (keep the original, just in case)
temp.list = strsplit(cdhitTE$class_family, split = "/")
testmp.list2 = lapply(temp.list, function(x) {as.data.frame(t(x))})
temp = do.call(rbind.fill, testmp.list2)
cdhitTE$Class = temp[,1]
cdhitTE$Family = temp[,2]
rm(temp.list, testmp.list2, temp)

lib.depth = sum(cdhitTE$depth)

#~~~~~~~~~~~~~~~~~~~~~~~~
# Assign Over-represented tags
#~~~~~~~~~~~~~~~~~~~~~~~~
# fit normal distribution to log transformed freq dist by maximum likelihood 
fd = fitdist(data = log10(cdhitTE$depth), distr = "norm")
# calculate outlier threshold 1 (3rd quartile + 1.5 IQR) in raw value
outlier.th1 = 10^(fd$estimate[1] + 2.68*fd$estimate[2] )

cdhitTE$Over = ifelse(cdhitTE$depth >=outlier.th1, "over", "no")

#~~~~~~~~~~~~~~~~~~~~~~~~
# Assign type by repStart
#~~~~~~~~~~~~~~~~~~~~~~~~
cdhitTE$RE_position = ifelse(cdhitTE$repStart >5, "external",  ifelse(is.na(cdhitTE$repStart),"none","internal"))
cdhitTE$RE_position[is.na(cdhitTE$RE_position)] = "none"

# save  to .RAD_TE.summary
write_csv(cdhitTE, paste0("./", substr(args[1],1,nchar(args[1])-9), ".RAD_TE.summary1"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate summary ONLY for over-represented tags, scaled against total lib.depth
#~~~~~~~~~~~~~~~~~~~~~~~~~~
cdhitTE_Over = subset(cdhitTE, Over == "over")
External = subset(cdhitTE_Over, RE_position=="external")
Internal =  subset(cdhitTE_Over, RE_position=="internal")

depth.summary = data.frame(sample = substr(args[1],1,nchar(args[1])-9),
                           # For all types
                           lib.depth = lib.depth,
                           DNA = sum(subset(cdhitTE_Over, Class=="DNA")$depth)/ lib.depth,
                           RC = sum(subset(cdhitTE_Over, Class=="RC")$depth)/ lib.depth,
                           LTR = sum(subset(cdhitTE_Over, Class=="LTR")$depth)/ lib.depth,
                           LINE = sum(subset(cdhitTE_Over, Class=="LINE")$depth)/ lib.depth,
                           SINE = sum(subset(cdhitTE_Over, Class=="SINE"|Class=="SINE?")$depth)/ lib.depth,
                           RNA = sum(subset(cdhitTE_Over, Class=="rRNA"|Class=="snRNA"|Class=="tRNA")$depth)/ lib.depth,
                           UO = sum(subset(cdhitTE_Over, Class=="Unknown"|Class=="Other")$depth)/ lib.depth,
                           Sim = sum(subset(cdhitTE_Over, Class=="Simple_repeat")$depth)/ lib.depth,
                           Sat = sum(subset(cdhitTE_Over, Class=="Satellite")$depth)/ lib.depth,
                           LC = sum(subset(cdhitTE_Over, Class=="Low_complexity")$depth)/ lib.depth,
                           # Internal
                           DNA.int = sum(subset(Internal, Class=="DNA")$depth)/ lib.depth,
                           RC.int = sum(subset(Internal, Class=="RC")$depth)/ lib.depth,
                           LTR.int = sum(subset(Internal, Class=="LTR")$depth)/ lib.depth,
                           LINE.int = sum(subset(Internal, Class=="LINE")$depth)/ lib.depth,
                           SINE.int = sum(subset(Internal, Class=="SINE"|Class=="SINE?")$depth)/ lib.depth,
                           SS.int = sum(subset(Internal, Class=="Simple_repeat"|Class=="Satellite")$depth)/ lib.depth,
                           RNA.int = sum(subset(Internal, Class=="rRNA"|Class=="snRNA"|Class=="tRNA")$depth)/ lib.depth,
                           UO.int = sum(subset(Internal, Class=="Unknown"|Class=="Other")$depth)/ lib.depth,
                           # External
                           DNA.ext = sum(subset(External, Class=="DNA")$depth)/ lib.depth,
                           RC.ext = sum(subset(External, Class=="RC")$depth)/ lib.depth,
                           LTR.ext = sum(subset(External, Class=="LTR")$depth)/ lib.depth,
                           LINE.ext = sum(subset(External, Class=="LINE")$depth)/ lib.depth,
                           SINE.ext = sum(subset(External, Class=="SINE"|Class=="SINE?")$depth)/ lib.depth,
                           SS.ext = sum(subset(External, Class=="Simple_repeat"|Class=="Satellite")$depth)/ lib.depth,
                           RNA.ext = sum(subset(External, Class=="rRNA"|Class=="snRNA"|Class=="tRNA")$depth)/ lib.depth,
                           UO.ext = sum(subset(External, Class=="Unknown"|Class=="Other")$depth)/ lib.depth,
                           # Overall
                           TE = sum(subset(cdhitTE_Over, Class=="DNA"|Class=="RC"|Class=="LTR"|Class=="LINE"|Class=="SINE"|Class=="SINE?"|Class=="Unknown"|Class=="Other")$depth)/ lib.depth,
                           TE.int = sum(subset(Internal, Class=="DNA"|Class=="RC"|Class=="LTR"|Class=="LINE"|Class=="SINE"|Class=="SINE?"|Class=="Unknown"|Class=="Other")$depth)/ lib.depth,
                           TE.ext = sum(subset(External, Class=="DNA"|Class=="RC"|Class=="LTR"|Class=="LINE"|Class=="SINE"|Class=="SINE?"|Class=="Unknown"|Class=="Other")$depth)/ lib.depth,
                           NoTE = sum(subset(cdhitTE_Over, RE_position=="none")$depth) / lib.depth
)
print("done")

write_csv(depth.summary, paste0("./", substr(args[1],1,nchar(args[1])-9), ".RAD_TE.summary2"))

  

