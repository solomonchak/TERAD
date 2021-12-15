#########################################
# RAD_TE_summary
#########################################
# updated 202101003

args <- commandArgs(trailingOnly = TRUE)

print(paste0("Input file 1: " ,args[1]))
print(paste0("Input file 2: " ,args[2]))
print(paste0("Input file 3: " ,args[3]))

# for testing purpose
# setwd("/Users/solomon/Desktop/Juan 20211214")
# args = c("L09_1_500k.fasta.cd.summary2", "L09_1_500k.fasta.cd.out", "L09_1_500k.fasta.cd.masked.annot")
# update.packages(ask = FALSE)

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

# In some case cdhit is very big, so it'll be better to remove the ones with low depth first. 
# Essentialy, everything under the third quantile could be removed first.
# But it's important to keep track of the total depth and the array of depth for later.
cdhit.depth.array = cdhit$depth
cdhit.total.depth = sum(cdhit.depth.array)

cdhit = cdhit[!cdhit$depth<quantile(cdhit$depth)[4],]



#~~~~~~~~~~~~~~~~~~~~~~~~
print("Read RepeatMasker's .cd.out")
#~~~~~~~~~~~~~~~~~~~~~~~~
RMout = read_table(paste0("./",args[2]), skip=3, col_types = "____cnn__cc____c", col_names = c("TE.tag", "repStart", "repEnd", "repeat", "class_family", "remove"), progress = F)

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
Protein = read_table(paste0("./",args[3]), skip=3, col_types = "c__cnn_cc__", col_names = c("pValue","TE.tag", "repStart", "repEnd", "repeat", "class_family"), progress = F)

# RMasker/TRF doesn't give a p-value
# They are being kep in the next step.

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

# cdhit is the largest dataframe and may be too large for R to handle. 
# Solution: break it up into smaller chuncks and process it serially to conserve memory
# 1000 rows seems to be a good amount to work with. Later version may allow user to specify this

chunk <- 1000
# chunk <- 10
n <- nrow(cdhit)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
cdhit_split <- split(cdhit,r)
# Test with just two lists
# cdhit_split = cdhit_split[c(1,2)]

# Combine RMout and Protein & remove duplicates
# there're some duplicates here that needs to be removed.
RMout_comb = ddply(rbind(RMout, Protein), "TE.tag18", function(x){
  x[1,]
})


# Do these analyses in each list of cdhit_split
fun1 = function(cdhit_element){
  # cdhit_element = cdhit_split[["1"]]
  
  cdhitTE = merge(cdhit_element, RMout_comb, by.x="rep18", by.y="TE.tag18", all = T)
  
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
  # cdhitTE$Family = temp[,2]
  cdhitTE$Family = tryCatch(temp[,2], error=function(err) NA) # In case all class_family are "Simple_repeat", then there'll be no second column.
  rm(temp.list, testmp.list2, temp)
  
  #lib.depth = sum(cdhitTE$depth)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~
  # Assign Over-represented tags ---- this needs to be done outside the function
  #~~~~~~~~~~~~~~~~~~~~~~~~
  # fit normal distribution to log transformed freq dist by maximum likelihood 
 # fd = fitdist(data = log10(cdhitTE$depth), distr = "norm")
  # calculate outlier threshold 1 (3rd quartile + 1.5 IQR) in raw value
  #outlier.th1 = 10^(fd$estimate[1] + 2.68*fd$estimate[2] )
  
 # cdhitTE$Over = ifelse(cdhitTE$depth >=outlier.th1, "over", "no")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~
  # Assign type by repStart
  #~~~~~~~~~~~~~~~~~~~~~~~~
  cdhitTE$RE_position = ifelse(cdhitTE$repStart >5, "external",  ifelse(is.na(cdhitTE$repStart),"none","internal"))
  cdhitTE$RE_position[is.na(cdhitTE$RE_position)] = "none"
  
  return(cdhitTE)
}
cdhitTE_list = lapply(cdhit_split, fun1)


#~~~~~~~~~~~~~~~~~~~~~~~~
# Assign Over-represented tags 
#~~~~~~~~~~~~~~~~~~~~~~~~
# Extract "depth" from the list
depth.array = unlist(lapply(cdhitTE_list, function(x){x$depth}))
# add the saved ones which were removed from cdhit
depth.array = c(depth.array, cdhit.depth.array)


# calculate total depth across the list
lib.depth = sum(depth.array)
# add the saved ones which were removed from cdhit
lib.depth = lib.depth + cdhit.total.depth

# calculate the threshold for assigning "Over"
# fit a normal distribution to log transformed freq dist by maximum likelihood 
fd = fitdist(data = log10(as.vector(depth.array)), distr = "norm")
# calculate outlier threshold 1 (3rd quartile + 1.5 IQR) in raw value
outlier.th1 = 10^(fd$estimate[1] + 2.68*fd$estimate[2] )


# Assign "Over" within list
cdhitTE_list = lapply(cdhitTE_list, function(x){
  x$Over = ifelse(x$depth >=outlier.th1, "over", "no")
  x
})


# save  to .RAD_TE.summary1.X
# Not really needed. It's here just in case it's needed somehow.
# sapply(names(cdhitTE_list), function (x) write_csv(cdhitTE_list[[x]], file=paste0("./", substr(args[1],1,nchar(args[1])-9), ".RAD_TE.summary1.", x) )   )




#~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate summary ONLY for over-represented tags, scaled against total lib.depth
#~~~~~~~~~~~~~~~~~~~~~~~~~~

# Trim the list to just "over" & then combine
cdhitTE_Over = do.call(rbind, lapply(cdhitTE_list, function(x){subset(x, Over =="over")}))
# there were no over-represented one in the test. Do more ###########################*******************************cont.
# the trick to reduce size is to remove all the under-represented one at the beginning of the script, then combine the df.


#cdhitTE_Over = subset(cdhitTE, Over == "over")
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

  

