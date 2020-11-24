# load packages
library(dplyr)
library(data.table)

# list bed files in out folder
variants.list = list.files(path = "out/bed/", 
                           recursive = TRUE,
                           pattern = "\\.bed$", 
                           full.names = TRUE)

#rea  d variants into R
variants = rbindlist(sapply(variants.list, 
                            fread, 
                            simplify = FALSE),
                     use.names = TRUE)

# compute variant size
variant.size = variants$V3 - variants$V2

# compute 2x of variant size
variant.size.x2 =  variant.size * 2

# compute new start and end position (adjascent regions flanking variants)
start = variants$V2 - variant.size.x2
end = variants$V3 + variant.size.x2

# construct new data frame
reads = variants

# input computed values
reads$V2 = start
reads$V3 = end

# compile new list
reads = reads %>%
  select(V1, V2, V3) %>% 
  mutate_if(is.numeric, round) %>%
  mutate(V2 = if_else(V2 < 0, 1, V2))

# subset on chr
chr1 = (subset(reads, V1 == "chr1"))
chr2 = (subset(reads, V1 == "chr2"))
chr3 = (subset(reads, V1 == "chr3"))
chr4 = (subset(reads, V1 == "chr4"))
chr5 = (subset(reads, V1 == "chr5"))
chr6 = (subset(reads, V1 == "chr6"))
chr7 = (subset(reads, V1 == "chr7"))
chr8 = (subset(reads, V1 == "chr8"))
chr9 = (subset(reads, V1 == "chr9"))
chr10 = (subset(reads, V1 == "chr10"))
chr11 = (subset(reads, V1 == "chr11"))
chr12 = (subset(reads, V1 == "chr12"))
chr13 = (subset(reads, V1 == "chr13"))
chr14 = (subset(reads, V1 == "chr14"))
chr15 = (subset(reads, V1 == "chr15"))
chr16 = (subset(reads, V1 == "chr16"))
chr17 = (subset(reads, V1 == "chr17"))
chr18 = (subset(reads, V1 == "chr18"))
chr19 = (subset(reads, V1 == "chr19"))
chr20 = (subset(reads, V1 == "chr20"))
chr21 = (subset(reads, V1 == "chr21"))
chr22 = (subset(reads, V1 == "chr22"))
chrx = (subset(reads, V1 == "chrx"))

# check if values are greater than chromosome end (GRCh38)
chr1 =  chr1 %>%
  mutate(V3 = ifelse(V3 > 248956422, 248956421, V3))

chr2 =  chr2 %>%
  mutate(V3 = ifelse(V3 > 242193529, 242193528, V3))

chr3 =  chr3 %>%
  mutate(V3 = ifelse(V3 > 198295559, 198295558, V3))

chr4 =  chr4 %>%
  mutate(V3 = ifelse(V3 > 190214555, 190214554, V3))

chr5 =  chr5 %>%
  mutate(V3 = ifelse(V3 > 181538259, 181538258, V3))

chr6 =  chr6 %>%
  mutate(V3 = ifelse(V3 > 170805979, 170805978, V3))

chr7 =  chr7 %>%
  mutate(V3 = ifelse(V3 > 159345973, 159345972, V3))

chr8 =  chr8 %>%
  mutate(V3 = ifelse(V3 > 145138636, 145138635, V3))

chr9 =  chr9 %>%
  mutate(V3 = ifelse(V3 > 138394717, 138394716, V3))

chr10 =  chr10 %>%
  mutate(V3 = ifelse(V3 > 133797422, 133797421, V3))

chr11 =  chr11 %>%
  mutate(V3 = ifelse(V3 > 135086622, 135086621, V3))

chr12 =  chr12 %>%
  mutate(V3 = ifelse(V3 > 133275309, 133275308, V3))

chr13 =  chr13 %>%
  mutate(V3 = ifelse(V3 > 114364328, 114364327, V3))

chr14 =  chr14 %>%
  mutate(V3 = ifelse(V3 > 107043718, 107043717, V3))

chr15 =  chr15 %>%
  mutate(V3 = ifelse(V3 > 101991189,101991188, V3))

chr16 =  chr16 %>%
  mutate(V3 = ifelse(V3 > 90338345, 90338344, V3))

chr17 =  chr17 %>%
  mutate(V3 = ifelse(V3 > 83257441, 83257440, V3))

chr18 =  chr18 %>%
  mutate(V3 = ifelse(V3 > 80373285, 80373284, V3))

chr19 =  chr19 %>%
  mutate(V3 = ifelse(V3 > 58617616, 58617615, V3))

chr20 =  chr20 %>%
  mutate(V3 = ifelse(V3 > 64444167, 64444166, V3))

chr21 =  chr21 %>%
  mutate(V3 = ifelse(V3 > 46709983, 100, V3))

chr22 =  chr22 %>%
  mutate(V3 = ifelse(V3 > 50818468, 50818467, V3))

chrx =  chrx %>%
  mutate(V3 = ifelse(V3 > 156040895, 156040894, V3))

# combine list
reads = rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrx)

# export as bed file
write.table(reads, "in/bed_reads/2x.adjascent.bed", quote = F, row.names = F, col.names = F, sep = "\t")