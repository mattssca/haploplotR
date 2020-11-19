# set sample name
sample = "HG00512"


# load packages
library(dplyr)

# read bedfiles
bedfile1 = read.table("../mattsada/Desktop/HGSV/NA192/NA19240/tables/H1_NA19240.txt", header = TRUE, sep = "\t")
bedfile2 = read.table("../mattsada/Desktop/HGSV/NA192/NA19240/tables/H2_NA19240.txt", header = TRUE, sep = "\t")

# convert back negative values for H1
bedfile1$start = bedfile1$start * -1
bedfile1$end = bedfile1$end * -1

# select variables for bedfile
bedfile1 = bedfile1 %>%
  select(chr, start, end, ID, probability, genotype)

# transform genotype to strand state  
levels(bedfile1$genotype)[levels(bedfile1$genotype)=="1|0"] = "+"
levels(bedfile1$genotype)[levels(bedfile1$genotype)=="1|1"] = "+"

# select variables for bedfile        
bedfile2 = bedfile2 %>%
  select(chr, start, end, ID, probability, genotype)

#transform genotype to strand state    
levels(bedfile2$genotype)[levels(bedfile2$genotype)=="0|1"] = "-"
levels(bedfile2$genotype)[levels(bedfile2$genotype)=="1|1"] = "-"
  
# bind both haplotypes into one bed file
bedfile = rbind(bedfile1, bedfile2)

# sort bedfile on chr
bedfile = bedfile[order(bedfile$chr),]
  
# export bedfile
write.table(bedfile, "out/bed/", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)