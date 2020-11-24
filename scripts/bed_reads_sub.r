# load packages
library(dplyr)
library(data.table)

# list bed files in out folder
variants.list = list.files(path = "out/bed/", 
                           recursive = TRUE,
                           pattern = "\\.bed$", 
                           full.names = TRUE)
#read variants into R
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


# export as bed file
write.table(reads, "in/bed_reads/2x.adjascent.bed", quote = F, row.names = F, col.names = F, sep = "\t")
