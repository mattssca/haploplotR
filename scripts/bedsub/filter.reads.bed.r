library(dplyr)

##########################################################HG005#####################################################################
# read variants intop R
HG00512.variants = read.table("Desktop/NEW BED/VARIANTS/hg00512.txt", sep = "\t", header = F)

# compute 50% of variant size
HG00512.50 =  HG00512.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
HG00512.start = HG00512.variants$V2 - HG00512.50
HG00512.end = HG00512.variants$V3 + HG00512.50

# construct new data frame
HG00512.reads = HG00512.variants

# input computed values
HG00512.reads$V2 = HG00512.start
HG00512.reads$V3 = HG00512.end

# compile new list
HG00512.reads = HG00512.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
  write.table(HG00512.reads, "Desktop/NEW BED/adjascent.regions/adj.hg00512.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# read variants intop R
HG00513.variants = read.table("Desktop/NEW BED/VARIANTS/hg00513.txt", sep = "\t", header = F)

# compute 50% of variant size
HG00513.50 =  HG00513.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
HG00513.start = HG00513.variants$V2 - HG00513.50
HG00513.end = HG00513.variants$V3 + HG00513.50

# construct new data frame
HG00513.reads = HG00513.variants

# input computed values
HG00513.reads$V2 = HG00513.start
HG00513.reads$V3 = HG00513.end

# compile new list
HG00513.reads = HG00513.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
write.table(HG00513.reads, "Desktop/NEW BED/adjascent.regions/adj.hg00513.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# read variants intop R
HG00514.variants = read.table("Desktop/NEW BED/VARIANTS/hg00514.txt", sep = "\t", header = F)

# compute 50% of variant size
HG00514.50 =  HG00514.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
HG00514.start = HG00514.variants$V2 - HG00514.50
HG00514.end = HG00514.variants$V3 + HG00514.50

# construct new data frame
HG00514.reads = HG00514.variants

# input computed values
HG00514.reads$V2 = HG00514.start
HG00514.reads$V3 = HG00514.end

# compile new list
HG00514.reads = HG00514.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
write.table(HG00514.reads, "Desktop/NEW BED/adjascent.regions/adj.hg00514.bed", quote = F, row.names = F, col.names = F, sep = "\t")

##########################################################HG007#####################################################################
# read variants intop R
HG00731.variants = read.table("Desktop/NEW BED/VARIANTS/hg00731.txt", sep = "\t", header = F)

# compute 50% of variant size
HG00731.50 =  HG00731.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
HG00731.start = HG00731.variants$V2 - HG00731.50
HG00731.end = HG00731.variants$V3 + HG00731.50

# construct new data frame
HG00731.reads = HG00731.variants

# input computed values
HG00731.reads$V2 = HG00731.start
HG00731.reads$V3 = HG00731.end

# compile new list
HG00731.reads = HG00731.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
  write.table(HG00731.reads, "Desktop/NEW BED/adjascent.regions/adj.HG00731.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# read variants intop R
HG00732.variants = read.table("Desktop/NEW BED/VARIANTS/hg00732.txt", sep = "\t", header = F)

# compute 50% of variant size
HG00732.50 =  HG00732.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
HG00732.start = HG00732.variants$V2 - HG00732.50
HG00732.end = HG00732.variants$V3 + HG00732.50

# construct new data frame
HG00732.reads = HG00732.variants

# input computed values
HG00732.reads$V2 = HG00732.start
HG00732.reads$V3 = HG00732.end

# compile new list
HG00732.reads = HG00732.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)


# export as bed file
write.table(HG00732.reads, "Desktop/NEW BED/adjascent.regions/adj.HG00732.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# read variants intop R
HG00733.variants = read.table("Desktop/NEW BED/VARIANTS/hg00733.txt", sep = "\t", header = F)

# compute 50% of variant size
HG00733.50 =  HG00733.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
HG00733.start = HG00733.variants$V2 - HG00733.50
HG00733.end = HG00733.variants$V3 + HG00733.50

# construct new data frame
HG00733.reads = HG00733.variants

# input computed values
HG00733.reads$V2 = HG00733.start
HG00733.reads$V3 = HG00733.end

# compile new list
HG00733.reads = HG00733.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
write.table(HG00733.reads, "Desktop/NEW BED/adjascent.regions/adj.HG00733.bed", quote = F, row.names = F, col.names = F, sep = "\t")

##########################################################NA192#####################################################################
# read variants intop R
NA19238.variants = read.table("Desktop/NEW BED/VARIANTS/na19238.txt", sep = "\t", header = F)

# compute 50% of variant size
NA19238.50 =  NA19238.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
NA19238.start = NA19238.variants$V2 - NA19238.50
NA19238.end = NA19238.variants$V3 + NA19238.50

# construct new data frame
NA19238.reads = NA19238.variants

# input computed values
NA19238.reads$V2 = NA19238.start
NA19238.reads$V3 = NA19238.end

# compile new list
NA19238.reads = NA19238.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
  write.table(NA19238.reads, "Desktop/NEW BED/adjascent.regions/adj.NA19238.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# read variants intop R
NA19239.variants = read.table("Desktop/NEW BED/VARIANTS/na19239.txt", sep = "\t", header = F)

# compute 50% of variant size
NA19239.50 =  NA19239.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
NA19239.start = NA19239.variants$V2 - NA19239.50
NA19239.end = NA19239.variants$V3 + NA19239.50

# construct new data frame
NA19239.reads = NA19239.variants

# input computed values
NA19239.reads$V2 = NA19239.start
NA19239.reads$V3 = NA19239.end

# compile new list
NA19239.reads = NA19239.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
write.table(NA19239.reads, "Desktop/NEW BED/adjascent.regions/adj.NA19239.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# read variants intop R
NA19240.variants = read.table("Desktop/NEW BED/VARIANTS/na19240.txt", sep = "\t", header = F)

# compute 50% of variant size
NA19240.50 =  NA19240.variants$V4 / 2

# compute new start and end position (adjascent regions flanking variants)
NA19240.start = NA19240.variants$V2 - NA19240.50
NA19240.end = NA19240.variants$V3 + NA19240.50

# construct new data frame
NA19240.reads = NA19240.variants

# input computed values
NA19240.reads$V2 = NA19240.start
NA19240.reads$V3 = NA19240.end

# compile new list
NA19240.reads = NA19240.reads %>%
  select(V1, V2, V3) %>% 
 mutate_if(is.numeric, round)

# export as bed file
write.table(NA19240.reads, "Desktop/NEW BED/adjascent.regions/adj.NA19240.bed", quote = F, row.names = F, col.names = F, sep = "\t")