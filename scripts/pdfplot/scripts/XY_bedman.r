bedfile1 = read.table("../mattsada/Desktop/HGSV/NA192/NA19239/tables/H1_NA19239.txt", header = TRUE, sep = "\t")
bedfile2 = read.table("../mattsada/Desktop/HGSV/NA192/NA19239/tables/H2_NA19239.txt", header = TRUE, sep = "\t")
bedfile3 = read.table("../mattsada/Desktop/HGSV/NA192/NA19239/tables/Haploid_NA19239.txt", header = TRUE, sep = "\t")

bedfile1$start = bedfile1$start * -1
bedfile1$end = bedfile1$end * -1

bedfile1 = bedfile1 %>%
  select(chr, start, end, ID, probability, genotype)

levels(bedfile1$genotype)[levels(bedfile1$genotype)=="1|0"] = "+"
levels(bedfile1$genotype)[levels(bedfile1$genotype)=="1|1"] = "+"

bedfile2 = bedfile2 %>%
  select(chr, start, end, ID, probability, genotype)

levels(bedfile2$genotype)[levels(bedfile2$genotype)=="0|1"] = "-"
levels(bedfile2$genotype)[levels(bedfile2$genotype)=="1|1"] = "-"

bedfile3 = bedfile3 %>%
  select(chr, start, end, ID, probability, genotype)

bedfile3$genotype = as.factor(bedfile3$genotype)
levels(bedfile3$genotype)[levels(bedfile3$genotype)=="1"] = "+"
levels(bedfile3$genotype)[levels(bedfile3$genotype)=="1"] = "+"

bedfile = rbind(bedfile1, bedfile2, bedfile3)

bedfile = bedfile[order(bedfile$chr),]

write.table(bedfile, "../mattsada/Desktop/pdfplot/in/BED/NA19239.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  