  # disable warnings
  options(warn=-1)
  
  # If a package is installed, it will be loaded. If any are not, the missing package(s) will be installed from CRAN and then loaded.
  # packages of interest
  packages = c("dplyr", "gridExtra", "ggplot2", "data.table", "psych")
  
  # load or install&load all
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )
  
  ################################################## define palette ########################################################
  dataset.color = ""
  homozygous.varaints.color = "black" 
  heterozygous.variants.color = "black"
  
  ################################################## data wrangling ########################################################
  
  # read data into R
  # List txt files
  variants.list = list.files(path = "./in/", 
                             recursive = TRUE,
                             pattern = "\\.txt$", 
                             full.names = TRUE)
  
  # Read all the files and create a ID column to store filenames
  variants = rbindlist(sapply(variants.list, 
                              fread, 
                              simplify = FALSE),
                       use.names = TRUE,
                       idcol = "ID" )
  
  # Keep only filenames (not full path in ID varaible)
  variants$ID <- gsub('./in//', '', variants$ID)
  variants$ID <- gsub('.txt', '', variants$ID)
  
  # set ID variable as factor
  variants$ID = as.factor(variants$ID)
  
  # set sample name
  sample.name = as.character(variants$ID[1], quote = FALSE)
  
  #subset on haplotype
  haplotype1 = variants %>% filter(genotype %in% c("1|0", "1|1"))
  haplotype2 = variants %>% filter(genotype %in% c("0|1", "1|1"))
  haploids = variants %>% filter(genotype %in% c("1"))
  
  # filter out events residing in black listed regions
  haplotype1 = haplotype1 %>% filter(!probability %in% c("blacklisted"))
  haplotype2 = haplotype2 %>% filter(!probability %in% c("blacklisted"))
  haploids = haploids %>% filter(!probability %in% c("blacklisted"))
  
  # filter out events with low probability score
  haplotype1 = subset(haplotype1, probability > 0.95)
  haplotype2 = subset(haplotype2, probability > 0.95)
  haploids = subset(haploids, probability > 0.95)
  
  # transform coordinates in haplotype1 to negative values
  haplotype1$start = haplotype1$start * -1
  haplotype1$end = haplotype1$end * -1
  
  # set variables for exporting
  out1 = sprintf("./out/tables/H1_%s.txt", sample.name)
  out2 = sprintf("./out/tables/H2_%s.txt", sample.name)
  outhap = sprintf("./out/tables/Haploid_%s.txt", sample.name)
  
  # export tables
  write.table(haplotype1, out1, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(haplotype2, out2, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(haploids, outhap, sep = "\t", quote = FALSE, row.names = FALSE)
  
  print("Sucess! Haplotypes exported to ./out")
  
  ################################################## plot ideograms ########################################################
   
  # read variants into R
  h1 = read.table(out1, sep = "\t", header = TRUE)
  h2 = read.table(out2, sep = "\t", header = TRUE)
  h = read.table(outhap, sep = "\t", header = TRUE)
  
  # subset both ahplotypes on genotype
  hetero.h1 = h1 %>% filter(genotype %in% c("1|0"))
  homo.h1 = h1 %>% filter(genotype %in% c("1|1"))
  hetero.h2 = h2 %>% filter(genotype %in% c("0|1"))
  homo.h2 = h2 %>% filter(genotype %in% c("1|1"))
  
  ################################################## generate bed files ########################################################
  
  # compile genome-wide bed files
  hbed1 = haplotype1 %>%
    select(chr, start, end, ID, probability, genotype)
  
  hbed2 = haplotype2 %>%
    select(chr, start, end, ID, probability, genotype)
  
  hbed = haploids %>%
    select(chr, start, end, ID, probability, genotype)
  
  # convert back negative values for haplotype1
  hbed1$start = hbed1$start * -1
  hbed1$end = hbed1$end * -1
  
  # convert genotype to factor
  hbed1$genotype = as.factor(hbed1$genotype)
  hbed2$genotype = as.factor(hbed2$genotype)
  hbed$genotype = as.factor(hbed$genotype)
  
  # convert chr to factor
  hbed1$chr = as.factor(hbed1$chr)
  hbed2$chr = as.factor(hbed2$chr)
  hbed$chr = as.factor(hbed$chr)
  
  # adjust genotype for bedfile
  levels(hbed1$genotype)[levels(hbed1$genotype)=="1|0"] = "+"
  levels(hbed1$genotype)[levels(hbed1$genotype)=="1|1"] = "+"
  levels(hbed2$genotype)[levels(hbed2$genotype)=="0|1"] = "-"
  levels(hbed2$genotype)[levels(hbed2$genotype)=="1|1"] = "-"
  levels(hbed$genotype)[levels(hbed$genotype)=="1"] = "-"
  
  # rbind both haplotypes into one bed file
  bedfile = rbind(hbed1, hbed2, hbed)
  
  # select chr, start, end variables for bed file fomrmat
  # sort bed file
  bedfile = bedfile %>% 
    select(chr, start, end, ID, probability, genotype)
  
  # calcualte variant size
  bedfile = bedfile %>% 
    mutate(size = end-start)
  
  # sort bed file
  bedfile = bedfile %>% 
    arrange(chr, start, end)
  
  # create summary tables
  bed.size = describeBy(bedfile$size)
  
  bed.size = bed.size %>% 
    select(n, mean, sd, median, trimmed, mad, min, max, range, skew, kurtosis, se)
  
  bedfile.out = bedfile %>% 
    select(chr, start, end, ID, probability, genotype)


  bedfile.out$start <- format(bedfile.out$start,scientific=FALSE)
  bedfile.out$end <- format(bedfile.out$end,scientific=FALSE)

  # set variables for exporting
  out4 = sprintf("./out/bed/%s.bed", sample.name)
  out5 = sprintf("./out/metrics/%s.summary.txt", sample.name)
  
  # export tables
  write.table(bedfile.out, out4, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(bed.size, out5, sep = "\t", quote = FALSE, row.names = FALSE)
  
  ################################################## variant size plotting ########################################################
  
  # add dataset variable to dataframe
  SV.calls = bedfile
  namevector <- c("cell.set")
  SV.calls[,namevector] <- sample.name
  
  #convert variable from numeric to factor
  SV.calls$cell.set = as.factor(SV.calls$cell.set)
  
  # plotting
  SV.size.violine = ggplot(SV.calls, aes(x = SV.calls$cell.set, y = SV.calls$size, fill = SV.calls$cell.set)) + 
    stat_summary(fun = mean, geom = "point", shape = 20, size = 3) +
    labs(title = "Inversion Size Distribution", subtitle = "GRCh38", x = "", y = "Size (bp)", fill = "Data") +
    theme(legend.position = "none") +
    scale_y_log10() +
    geom_violin(trim = FALSE, scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), fill = "grey46")
  
  # binned varaint sizes
  # set up cut-off values 
  breaks <- c(1, 50, 1000, 10000, 100000, 1000000, 5000000)
  
  # specify interval/bin labels
  tags <- c("1-50bp","50bp-1kb", "1kb-10kb", "10kb-100kb", "100kb-1Mb", "1Mb-5Mb")
  
  # bucketing values into bins
  SV_tags <- cut(SV.calls$size, breaks=breaks, include.lowest=TRUE, right=FALSE, labels=tags)
  
  # plotting
  SV.binned = ggplot(data = as_tibble(SV_tags), mapping = aes(x=value)) + 
    geom_bar(color = "black", fill = "grey46") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) +
    labs(title = "Binned Variant Size", subtitle = "GRCh38", x = "", y = "Size (bp)", fill = "")
  
  # variants per chromosome counts
  # add chr count variable
  SV.calls = SV.calls %>% 
    add_count(chr)
  
  # subset dataframes to only include chr and count (n) + call.set
  SV.calls.sub = select(SV.calls, chr, n, cell.set)
  
  # remove duplicate chr rows
  SV.calls.sub = distinct(SV.calls.sub, chr, .keep_all = TRUE)
  
  # set max limmits for y scale
  ymax = max(SV.calls.sub$n)
  
  # plotting
  SV.chrdist.box = ggplot(SV.calls.sub, aes(x = SV.calls.sub$chr, y = SV.calls.sub$n, fill = SV.calls.sub$cell.set)) +
    labs(title = "Variants per Chromosome", subtitle = "", x = "", y = "Variant count", fill = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(), axis.title.x = element_text(), axis.title.y = element_text(), plot.title = element_text(), plot.subtitle = element_text(), legend.title = element_text(), legend.text = element_text(), legend.position = "top") +
    scale_x_discrete(limits=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) +
    geom_bar(position = position_dodge(), stat = "identity", color = "grey3") +
    scale_y_continuous(limits = c(0, ymax), breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)) +
    scale_fill_manual(values=c("grey46"))
  
  # export plots
  ggsave("viol.png", SV.size.violine, path = "./out/figs/", limitsize = FALSE, width = 10, height = 10, units = c("cm"), scale = 1, dpi = 300)
  ggsave("binned.box.png", SV.binned, path = "./out/figs/", limitsize = FALSE, width = 10, height = 10, units = c("cm"), scale = 1, dpi = 300)
  ggsave("chrdist.png", SV.chrdist.box, path = "./out/figs/", limitsize = FALSE, width = 20, height = 10, units = c("cm"), scale = 1, dpi = 300)
  
  print("plots sucessfully generated!")
  
  ################################# ideograms ################################
  # Subset data and prepare for plotting
  ideo.hetero.h1 = h1 %>% filter(genotype %in% c("1|0"))
  ideo.homo.h1 = h1 %>% filter(genotype %in% c("1|1"))
  ideo.hetero.h2 = h2 %>% filter(genotype %in% c("0|1"))
  ideo.homo.h2 = h2 %>% filter(genotype %in% c("1|1"))
  ideo.haploid = h %>% filter(genotype %in% c("1"))
  
  # transform negative values for H1 coordinates
  ideo.hetero.h1$start = ideo.hetero.h1$start * -1
  ideo.hetero.h1$end = ideo.hetero.h1$end * -1
  ideo.homo.h1$start = ideo.homo.h1$start * -1
  ideo.homo.h1$end = ideo.homo.h1$end * -1
  
  # Select varaibles for new data frame
  ideo.hetero.h1 = ideo.hetero.h1 %>%
    select(chr, start, end, genotype)
  
  ideo.homo.h1 = ideo.homo.h1 %>%
    select(chr, start, end, genotype)
  
  ideo.hetero.h2 = ideo.hetero.h2 %>%
    select(chr, start, end, genotype)
  
  ideo.homo.h2 = ideo.homo.h2 %>%
    select(chr, start, end, genotype)
  
  ideo.haploid = ideo.haploid %>%
    select(chr, start, end, genotype)
  
  # Transform genotype values to numeric values for plotting along y axis (lollipops)
  levels(ideo.hetero.h1$genotype)[levels(ideo.hetero.h1$genotype)=="1|0"] = "3"
  levels(ideo.homo.h1$genotype)[levels(ideo.homo.h1$genotype)=="1|1"] = "3"
  levels(ideo.hetero.h2$genotype)[levels(ideo.hetero.h2$genotype)=="0|1"] = "-3"
  levels(ideo.homo.h2$genotype)[levels(ideo.homo.h2$genotype)=="1|1"] = "-3"
  ideo.haploid$genotype <- sub("1", "2", ideo.haploid$genotype)
  
  # convert factors to numeric values
  ideo.hetero.h1$genotype = as.numeric(as.character(ideo.hetero.h1$genotype))
  ideo.homo.h1$genotype = as.numeric(as.character(ideo.homo.h1$genotype))
  ideo.hetero.h2$genotype = as.numeric(as.character(ideo.hetero.h2$genotype))
  ideo.homo.h2$genotype = as.numeric(as.character(ideo.homo.h2$genotype))
  ideo.haploid$genotype = as.numeric(as.character(ideo.haploid$genotype))
  
  # Transform coordinates to be in middle of each inversion
  ideo.hetero.h1 = ideo.hetero.h1 %>%
    mutate(size = end - start) %>%
    mutate(middle = size/2) %>%
    mutate(coordinate = start + middle) %>%
    mutate(genotype = genotype) %>%
    select(chr, start, end, size, coordinate, genotype)
  
  ideo.hetero.h2 = ideo.hetero.h2 %>%
    mutate(size = end - start) %>%
    mutate(middle = size/2) %>%
    mutate(coordinate = start + middle) %>%
    mutate(genotype = genotype) %>%
    select(chr, start, end, size, coordinate, genotype)
  
  ideo.homo.h1 = ideo.homo.h1 %>%
    mutate(size = end - start) %>%
    mutate(middle = size/2) %>%
    mutate(coordinate = start + middle) %>%
    mutate(genotype = genotype) %>%
    select(chr, start, end, size, coordinate, genotype)
  
  ideo.homo.h2 = ideo.homo.h2 %>%
    mutate(size = end - start) %>%
    mutate(middle = size/2) %>%
    mutate(coordinate = start + middle) %>%
    mutate(genotype = genotype) %>%
    select(chr, start, end, size, coordinate, genotype)
  
  ideo.haploid = ideo.haploid %>%
    mutate(size = end - start) %>%
    mutate(middle = size/2) %>%
    mutate(coordinate = start + middle) %>%
    mutate(genotype = genotype) %>%
    select(chr, start, end, size, coordinate, genotype)
  
  # subset larg inversions (>10kb)
  ideo.hetero.h1.large = subset(ideo.hetero.h1, size > 50000)
  ideo.homo.h1.large = subset(ideo.homo.h1, size > 50000)
  ideo.hetero.h2.large = subset(ideo.hetero.h2, size > 50000)
  ideo.homo.h2.large = subset(ideo.homo.h2, size > 50000)
  ideo.haploid.large = subset(ideo.haploid, size > 50000)
  
  # subset data on chromosome
  # H1 Heterozygous
  chr1.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr1")
  chr2.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr2")
  chr3.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr3")
  chr4.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr4")
  chr5.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr5")
  chr6.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr6")
  chr7.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr7")
  chr8.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr8")
  chr9.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr9")
  chr10.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr10")
  chr11.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr11")
  chr12.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr12")
  chr13.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr13")
  chr14.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr14")
  chr15.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr15")
  chr16.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr16")
  chr17.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr17")
  chr18.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr18")
  chr19.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr19")
  chr20.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr20")
  chr21.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr21")
  chr22.ideo.h1.het = filter(ideo.hetero.h1, chr == "chr22")
  chrx.ideo = filter(ideo.haploid, chr == "chrX")
  chry.ideo = filter(ideo.haploid, chr == "chrY")
  
  # H1 heterozygous large
  chr1.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr1")
  chr2.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr2")
  chr3.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr3")
  chr4.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr4")
  chr5.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr5")
  chr6.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr6")
  chr7.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr7")
  chr8.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr8")
  chr9.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr9")
  chr10.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr10")
  chr11.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr11")
  chr12.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr12")
  chr13.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr13")
  chr14.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr14")
  chr15.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr15")
  chr16.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr16")
  chr17.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr17")
  chr18.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr18")
  chr19.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr19")
  chr20.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr20")
  chr21.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr21")
  chr22.ideo.h1.het.large = filter(ideo.hetero.h1.large, chr == "chr22")
  chrx.ideo.haploid.large = filter(ideo.haploid.large, chr == "chrX")
  chry.ideo.haploid.large = filter(ideo.haploid.large, chr == "chrY")
  
  # H1 Homozygous
  chr1.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr1")
  chr2.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr2")
  chr3.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr3")
  chr4.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr4")
  chr5.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr5")
  chr6.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr6")
  chr7.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr7")
  chr8.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr8")
  chr9.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr9")
  chr10.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr10")
  chr11.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr11")
  chr12.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr12")
  chr13.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr13")
  chr14.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr14")
  chr15.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr15")
  chr16.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr16")
  chr17.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr17")
  chr18.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr18")
  chr19.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr19")
  chr20.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr20")
  chr21.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr21")
  chr22.ideo.h1.hom = filter(ideo.homo.h1, chr == "chr22")

  # H1 homozygous large
  chr1.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr1")
  chr2.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr2")
  chr3.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr3")
  chr4.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr4")
  chr5.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr5")
  chr6.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr6")
  chr7.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr7")
  chr8.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr8")
  chr9.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr9")
  chr10.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr10")
  chr11.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr11")
  chr12.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr12")
  chr13.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr13")
  chr14.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr14")
  chr15.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr15")
  chr16.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr16")
  chr17.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr17")
  chr18.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr18")
  chr19.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr19")
  chr20.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr20")
  chr21.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr21")
  chr22.ideo.h1.hom.large = filter(ideo.homo.h1.large, chr == "chr22")

  # H2 Heterozygous
  chr1.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr1")
  chr2.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr2")
  chr3.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr3")
  chr4.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr4")
  chr5.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr5")
  chr6.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr6")
  chr7.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr7")
  chr8.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr8")
  chr9.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr9")
  chr10.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr10")
  chr11.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr11")
  chr12.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr12")
  chr13.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr13")
  chr14.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr14")
  chr15.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr15")
  chr16.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr16")
  chr17.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr17")
  chr18.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr18")
  chr19.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr19")
  chr20.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr20")
  chr21.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr21")
  chr22.ideo.h2.het = filter(ideo.hetero.h2, chr == "chr22")

  # H2 heterozygous large
  chr1.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr1")
  chr2.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr2")
  chr3.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr3")
  chr4.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr4")
  chr5.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr5")
  chr6.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr6")
  chr7.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr7")
  chr8.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr8")
  chr9.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr9")
  chr10.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr10")
  chr11.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr11")
  chr12.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr12")
  chr13.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr13")
  chr14.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr14")
  chr15.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr15")
  chr16.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr16")
  chr17.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr17")
  chr18.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr18")
  chr19.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr19")
  chr20.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr20")
  chr21.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr21")
  chr22.ideo.h2.het.large = filter(ideo.hetero.h2.large, chr == "chr22")

  # H2 Homozygous
  chr1.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr1")
  chr2.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr2")
  chr3.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr3")
  chr4.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr4")
  chr5.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr5")
  chr6.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr6")
  chr7.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr7")
  chr8.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr8")
  chr9.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr9")
  chr10.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr10")
  chr11.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr11")
  chr12.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr12")
  chr13.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr13")
  chr14.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr14")
  chr15.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr15")
  chr16.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr16")
  chr17.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr17")
  chr18.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr18")
  chr19.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr19")
  chr20.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr20")
  chr21.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr21")
  chr22.ideo.h2.hom = filter(ideo.homo.h2, chr == "chr22")

  # H2 homozygous large
  chr1.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr1")
  chr2.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr2")
  chr3.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr3")
  chr4.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr4")
  chr5.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr5")
  chr6.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr6")
  chr7.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr7")
  chr8.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr8")
  chr9.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr9")
  chr10.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr10")
  chr11.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr11")
  chr12.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr12")
  chr13.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr13")
  chr14.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr14")
  chr15.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr15")
  chr16.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr16")
  chr17.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr17")
  chr18.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr18")
  chr19.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr19")
  chr20.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr20")
  chr21.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr21")
  chr22.ideo.h2.hom.large = filter(ideo.homo.h2.large, chr == "chr22")

  # read dependencies for plotting
  chr1.table = read.table("dep/ideograms/chrtable/GRCh38.chr1.table.txt", header = FALSE, sep = "\t")
  chr2.table = read.table("dep/ideograms/chrtable/GRCh38.chr2.table.txt", header = FALSE, sep = "\t")
  chr3.table = read.table("dep/ideograms/chrtable/GRCh38.chr3.table.txt", header = FALSE, sep = "\t")
  chr4.table = read.table("dep/ideograms/chrtable/GRCh38.chr4.table.txt", header = FALSE, sep = "\t")
  chr5.table = read.table("dep/ideograms/chrtable/GRCh38.chr5.table.txt", header = FALSE, sep = "\t")
  chr6.table = read.table("dep/ideograms/chrtable/GRCh38.chr6.table.txt", header = FALSE, sep = "\t")
  chr7.table = read.table("dep/ideograms/chrtable/GRCh38.chr7.table.txt", header = FALSE, sep = "\t")
  chr8.table = read.table("dep/ideograms/chrtable/GRCh38.chr8.table.txt", header = FALSE, sep = "\t")
  chr9.table = read.table("dep/ideograms/chrtable/GRCh38.chr9.table.txt", header = FALSE, sep = "\t")
  chr10.table = read.table("dep/ideograms/chrtable/GRCh38.chr10.table.txt", header = FALSE, sep = "\t")
  chr11.table = read.table("dep/ideograms/chrtable/GRCh38.chr11.table.txt", header = FALSE, sep = "\t")
  chr12.table = read.table("dep/ideograms/chrtable/GRCh38.chr12.table.txt", header = FALSE, sep = "\t")
  chr13.table = read.table("dep/ideograms/chrtable/GRCh38.chr13.table.txt", header = FALSE, sep = "\t")
  chr14.table = read.table("dep/ideograms/chrtable/GRCh38.chr14.table.txt", header = FALSE, sep = "\t")
  chr15.table = read.table("dep/ideograms/chrtable/GRCh38.chr15.table.txt", header = FALSE, sep = "\t")
  chr16.table = read.table("dep/ideograms/chrtable/GRCh38.chr16.table.txt", header = FALSE, sep = "\t")
  chr17.table = read.table("dep/ideograms/chrtable/GRCh38.chr17.table.txt", header = FALSE, sep = "\t")
  chr18.table = read.table("dep/ideograms/chrtable/GRCh38.chr18.table.txt", header = FALSE, sep = "\t")
  chr19.table = read.table("dep/ideograms/chrtable/GRCh38.chr19.table.txt", header = FALSE, sep = "\t")
  chr20.table = read.table("dep/ideograms/chrtable/GRCh38.chr20.table.txt", header = FALSE, sep = "\t")
  chr21.table = read.table("dep/ideograms/chrtable/GRCh38.chr21.table.txt", header = FALSE, sep = "\t")
  chr22.table = read.table("dep/ideograms/chrtable/GRCh38.chr22.table.txt", header = FALSE, sep = "\t")
  chrx.table = read.table("dep/ideograms/chrtable/GRCh38.chrx.table.txt", header = FALSE, sep = "\t")
  chry.table = read.table("dep/ideograms/chrtable/GRCh38.chry.table.txt", header = FALSE, sep = "\t")
  
  chr1.cent = read.table("dep/ideograms/centromeres/GRCh38.chr1.centromeres.txt", header = FALSE, sep = "\t")
  chr2.cent = read.table("dep/ideograms/centromeres/GRCh38.chr2.centromeres.txt", header = FALSE, sep = "\t")
  chr3.cent = read.table("dep/ideograms/centromeres/GRCh38.chr3.centromeres.txt", header = FALSE, sep = "\t")
  chr4.cent = read.table("dep/ideograms/centromeres/GRCh38.chr4.centromeres.txt", header = FALSE, sep = "\t")
  chr5.cent = read.table("dep/ideograms/centromeres/GRCh38.chr5.centromeres.txt", header = FALSE, sep = "\t")
  chr6.cent = read.table("dep/ideograms/centromeres/GRCh38.chr6.centromeres.txt", header = FALSE, sep = "\t")
  chr7.cent = read.table("dep/ideograms/centromeres/GRCh38.chr7.centromeres.txt", header = FALSE, sep = "\t")
  chr8.cent = read.table("dep/ideograms/centromeres/GRCh38.chr8.centromeres.txt", header = FALSE, sep = "\t")
  chr9.cent = read.table("dep/ideograms/centromeres/GRCh38.chr9.centromeres.txt", header = FALSE, sep = "\t")
  chr10.cent = read.table("dep/ideograms/centromeres/GRCh38.chr10.centromeres.txt", header = FALSE, sep = "\t")
  chr11.cent = read.table("dep/ideograms/centromeres/GRCh38.chr11.centromeres.txt", header = FALSE, sep = "\t")
  chr12.cent = read.table("dep/ideograms/centromeres/GRCh38.chr12.centromeres.txt", header = FALSE, sep = "\t")
  chr13.cent = read.table("dep/ideograms/centromeres/GRCh38.chr13.centromeres.txt", header = FALSE, sep = "\t")
  chr14.cent = read.table("dep/ideograms/centromeres/GRCh38.chr14.centromeres.txt", header = FALSE, sep = "\t")
  chr15.cent = read.table("dep/ideograms/centromeres/GRCh38.chr15.centromeres.txt", header = FALSE, sep = "\t")
  chr16.cent = read.table("dep/ideograms/centromeres/GRCh38.chr16.centromeres.txt", header = FALSE, sep = "\t")
  chr17.cent = read.table("dep/ideograms/centromeres/GRCh38.chr17.centromeres.txt", header = FALSE, sep = "\t")
  chr18.cent = read.table("dep/ideograms/centromeres/GRCh38.chr18.centromeres.txt", header = FALSE, sep = "\t")
  chr19.cent = read.table("dep/ideograms/centromeres/GRCh38.chr19.centromeres.txt", header = FALSE, sep = "\t")
  chr20.cent = read.table("dep/ideograms/centromeres/GRCh38.chr20.centromeres.txt", header = FALSE, sep = "\t")
  chr21.cent = read.table("dep/ideograms/centromeres/GRCh38.chr21.centromeres.txt", header = FALSE, sep = "\t")
  chr22.cent = read.table("dep/ideograms/centromeres/GRCh38.chr22.centromeres.txt", header = FALSE, sep = "\t")
  chrx.cent = read.table("dep/ideograms/centromeres/GRCh38.chrx.centromeres.txt", header = FALSE, sep = "\t")
  chry.cent = read.table("dep/ideograms/centromeres/GRCh38.chry.centromeres.txt", header = FALSE, sep = "\t")
  
  # calculate number of inversion per haplotype
  chr1.h1.n = nrow(chr1.ideo.h1.het) + nrow(chr1.ideo.h1.hom)
  chr1.h2.n = nrow(chr1.ideo.h2.het) + nrow(chr1.ideo.h2.hom)
  
  chr2.h1.n = nrow(chr2.ideo.h1.het) + nrow(chr2.ideo.h1.hom)
  chr2.h2.n = nrow(chr2.ideo.h2.het) + nrow(chr2.ideo.h2.hom)
  
  chr3.h1.n = nrow(chr3.ideo.h1.het) + nrow(chr3.ideo.h1.hom)
  chr3.h2.n = nrow(chr3.ideo.h2.het) + nrow(chr3.ideo.h2.hom)
  
  chr4.h1.n = nrow(chr4.ideo.h1.het) + nrow(chr4.ideo.h1.hom)
  chr4.h2.n = nrow(chr4.ideo.h2.het) + nrow(chr4.ideo.h2.hom)
  
  chr5.h1.n = nrow(chr5.ideo.h1.het) + nrow(chr5.ideo.h1.hom)
  chr5.h2.n = nrow(chr5.ideo.h2.het) + nrow(chr5.ideo.h2.hom)
  
  chr6.h1.n = nrow(chr6.ideo.h1.het) + nrow(chr6.ideo.h1.hom)
  chr6.h2.n = nrow(chr6.ideo.h2.het) + nrow(chr6.ideo.h2.hom)
  
  chr7.h1.n = nrow(chr7.ideo.h1.het) + nrow(chr7.ideo.h1.hom)
  chr7.h2.n = nrow(chr7.ideo.h2.het) + nrow(chr7.ideo.h2.hom)
  
  chr8.h1.n = nrow(chr8.ideo.h1.het) + nrow(chr8.ideo.h1.hom)
  chr8.h2.n = nrow(chr8.ideo.h2.het) + nrow(chr8.ideo.h2.hom)
  
  chr9.h1.n = nrow(chr9.ideo.h1.het) + nrow(chr9.ideo.h1.hom)
  chr9.h2.n = nrow(chr9.ideo.h2.het) + nrow(chr9.ideo.h2.hom)
  
  chr10.h1.n = nrow(chr10.ideo.h1.het) + nrow(chr10.ideo.h1.hom)
  chr10.h2.n = nrow(chr10.ideo.h2.het) + nrow(chr10.ideo.h2.hom)
  
  chr11.h1.n = nrow(chr11.ideo.h1.het) + nrow(chr11.ideo.h1.hom)
  chr11.h2.n = nrow(chr11.ideo.h2.het) + nrow(chr11.ideo.h2.hom)
  
  chr12.h1.n = nrow(chr12.ideo.h1.het) + nrow(chr12.ideo.h1.hom)
  chr12.h2.n = nrow(chr12.ideo.h2.het) + nrow(chr12.ideo.h2.hom)
  
  chr13.h1.n = nrow(chr13.ideo.h1.het) + nrow(chr13.ideo.h1.hom)
  chr13.h2.n = nrow(chr13.ideo.h2.het) + nrow(chr13.ideo.h2.hom)
  
  chr14.h1.n = nrow(chr14.ideo.h1.het) + nrow(chr14.ideo.h1.hom)
  chr14.h2.n = nrow(chr14.ideo.h2.het) + nrow(chr14.ideo.h2.hom)
  
  chr15.h1.n = nrow(chr15.ideo.h1.het) + nrow(chr15.ideo.h1.hom)
  chr15.h2.n = nrow(chr15.ideo.h2.het) + nrow(chr15.ideo.h2.hom)
  
  chr16.h1.n = nrow(chr16.ideo.h1.het) + nrow(chr16.ideo.h1.hom)
  chr16.h2.n = nrow(chr16.ideo.h2.het) + nrow(chr16.ideo.h2.hom)
  
  chr17.h1.n = nrow(chr17.ideo.h1.het) + nrow(chr17.ideo.h1.hom)
  chr17.h2.n = nrow(chr17.ideo.h2.het) + nrow(chr17.ideo.h2.hom)
  
  chr18.h1.n = nrow(chr18.ideo.h1.het) + nrow(chr18.ideo.h1.hom)
  chr18.h2.n = nrow(chr18.ideo.h2.het) + nrow(chr18.ideo.h2.hom)
  
  chr19.h1.n = nrow(chr19.ideo.h1.het) + nrow(chr19.ideo.h1.hom)
  chr19.h2.n = nrow(chr19.ideo.h2.het) + nrow(chr19.ideo.h2.hom)
  
  chr20.h1.n = nrow(chr20.ideo.h1.het) + nrow(chr20.ideo.h1.hom)
  chr20.h2.n = nrow(chr20.ideo.h2.het) + nrow(chr20.ideo.h2.hom)
  
  chr21.h1.n = nrow(chr21.ideo.h1.het) + nrow(chr21.ideo.h1.hom)
  chr21.h2.n = nrow(chr21.ideo.h2.het) + nrow(chr21.ideo.h2.hom)
  
  chr22.h1.n = nrow(chr22.ideo.h1.het) + nrow(chr22.ideo.h1.hom)
  chr22.h2.n = nrow(chr22.ideo.h2.het) + nrow(chr22.ideo.h2.hom)
  
  chrx.n = nrow(chrx.ideo)
  chry.n = nrow(chry.ideo)
  
  # check for na observation and replace with 0
  no.obs = data.frame(0, 0, 0, 0, -1000000000, 0)
  names(no.obs) <- c("chr", "start", "end", "size", "coordinate", "genotype")  
  
  if(nrow(chr1.ideo.h1.het) == 0){
    chr1.ideo.h1.het <- rbind(chr1.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr2.ideo.h1.het) == 0){
    chr2.ideo.h1.het <- rbind(chr2.ideo.h1.het, no.obs) 
  }

  if(nrow(chr3.ideo.h1.het) == 0){
    chr3.ideo.h1.het <- rbind(chr3.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr4.ideo.h1.het) == 0){
    chr4.ideo.h1.het <- rbind(chr4.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr5.ideo.h1.het) == 0){
    chr5.ideo.h1.het <- rbind(chr5.ideo.h1.het, no.obs) 
  }

  if(nrow(chr6.ideo.h1.het) == 0){
    chr6.ideo.h1.het <- rbind(chr6.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr7.ideo.h1.het) == 0){
    chr7.ideo.h1.het <- rbind(chr7.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr8.ideo.h1.het) == 0){
    chr8.ideo.h1.het <- rbind(chr8.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr9.ideo.h1.het) == 0){
    chr9.ideo.h1.het <- rbind(chr9.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr10.ideo.h1.het) == 0){
    chr10.ideo.h1.het <- rbind(chr10.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr11.ideo.h1.het) == 0){
    chr11.ideo.h1.het <- rbind(chr11.ideo.h1.het, no.obs) 
  }

  if(nrow(chr12.ideo.h1.het) == 0){
    chr12.ideo.h1.het <- rbind(chr12.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr13.ideo.h1.het) == 0){
    chr13.ideo.h1.het <- rbind(chr13.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr14.ideo.h1.het) == 0){
    chr14.ideo.h1.het <- rbind(chr14.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr15.ideo.h1.het) == 0){
    chr15.ideo.h1.het <- rbind(chr15.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr16.ideo.h1.het) == 0){
    chr16.ideo.h1.het <- rbind(chr16.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr17.ideo.h1.het) == 0){
    chr17.ideo.h1.het <- rbind(chr17.ideo.h1.het, no.obs) 
  }

  if(nrow(chr18.ideo.h1.het) == 0){
    chr18.ideo.h1.het <- rbind(chr18.ideo.h1.het, no.obs) 
  }

  if(nrow(chr19.ideo.h1.het) == 0){
    chr19.ideo.h1.het <- rbind(chr19.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr20.ideo.h1.het) == 0){
    chr20.ideo.h1.het <- rbind(chr20.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr21.ideo.h1.het) == 0){
    chr21.ideo.h1.het <- rbind(chr21.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chr22.ideo.h1.het) == 0){
    chr22.ideo.h1.het <- rbind(chr22.ideo.h1.het, no.obs) 
  }
  
  if(nrow(chrx.ideo) == 0){
    chrx.ideo <- rbind(chrx.ideo, no.obs) 
  }
  
  if(nrow(chry.ideo) == 0){
    chry.ideo <- rbind(chry.ideo, no.obs) 
  }
  
  if(nrow(chr1.ideo.h1.hom) == 0){
    chr1.ideo.h1.hom <- rbind(chr1.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr2.ideo.h1.hom) == 0){
    chr2.ideo.h1.hom <- rbind(chr2.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr3.ideo.h1.hom) == 0){
    chr3.ideo.h1.hom <- rbind(chr3.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr4.ideo.h1.hom) == 0){
    chr4.ideo.h1.hom <- rbind(chr4.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr5.ideo.h1.hom) == 0){
    chr5.ideo.h1.hom <- rbind(chr5.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr6.ideo.h1.hom) == 0){
    chr6.ideo.h1.hom <- rbind(chr6.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr7.ideo.h1.hom) == 0){
    chr7.ideo.h1.hom <- rbind(chr7.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr8.ideo.h1.hom) == 0){
    chr8.ideo.h1.hom <- rbind(chr8.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr9.ideo.h1.hom) == 0){
    chr9.ideo.h1.hom <- rbind(chr9.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr10.ideo.h1.hom) == 0){
    chr10.ideo.h1.hom <- rbind(chr10.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr11.ideo.h1.hom) == 0){
    chr11.ideo.h1.hom <- rbind(chr11.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr12.ideo.h1.hom) == 0){
    chr12.ideo.h1.hom <- rbind(chr12.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr13.ideo.h1.hom) == 0){
    chr13.ideo.h1.hom <- rbind(chr13.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr14.ideo.h1.hom) == 0){
    chr14.ideo.h1.hom <- rbind(chr14.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr15.ideo.h1.hom) == 0){
    chr15.ideo.h1.hom <- rbind(chr15.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr16.ideo.h1.hom) == 0){
    chr16.ideo.h1.hom <- rbind(chr16.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr17.ideo.h1.hom) == 0){
    chr17.ideo.h1.hom <- rbind(chr17.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr18.ideo.h1.hom) == 0){
    chr18.ideo.h1.hom <- rbind(chr18.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr19.ideo.h1.hom) == 0){
    chr19.ideo.h1.hom <- rbind(chr19.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr20.ideo.h1.hom) == 0){
    chr20.ideo.h1.hom <- rbind(chr20.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr21.ideo.h1.hom) == 0){
    chr21.ideo.h1.hom <- rbind(chr21.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr22.ideo.h1.hom) == 0){
    chr22.ideo.h1.hom <- rbind(chr22.ideo.h1.hom, no.obs) 
  }
  
  if(nrow(chr1.ideo.h2.het) == 0){
    chr1.ideo.h2.het <- rbind(chr1.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr2.ideo.h2.het) == 0){
    chr2.ideo.h2.het <- rbind(chr2.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr3.ideo.h2.het) == 0){
    chr3.ideo.h2.het <- rbind(chr3.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr4.ideo.h2.het) == 0){
    chr4.ideo.h2.het <- rbind(chr4.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr5.ideo.h2.het) == 0){
    chr5.ideo.h2.het <- rbind(chr5.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr6.ideo.h2.het) == 0){
    chr6.ideo.h2.het <- rbind(chr6.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr7.ideo.h2.het) == 0){
    chr7.ideo.h2.het <- rbind(chr7.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr8.ideo.h2.het) == 0){
    chr8.ideo.h2.het <- rbind(chr8.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr9.ideo.h2.het) == 0){
    chr9.ideo.h2.het <- rbind(chr9.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr10.ideo.h2.het) == 0){
    chr10.ideo.h2.het <- rbind(chr10.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr11.ideo.h2.het) == 0){
    chr11.ideo.h2.het <- rbind(chr11.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr12.ideo.h2.het) == 0){
    chr12.ideo.h2.het <- rbind(chr12.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr13.ideo.h2.het) == 0){
    chr13.ideo.h2.het <- rbind(chr13.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr14.ideo.h2.het) == 0){
    chr14.ideo.h2.het <- rbind(chr14.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr15.ideo.h2.het) == 0){
    chr15.ideo.h2.het <- rbind(chr15.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr16.ideo.h2.het) == 0){
    chr16.ideo.h2.het <- rbind(chr16.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr17.ideo.h2.het) == 0){
    chr17.ideo.h2.het <- rbind(chr17.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr18.ideo.h2.het) == 0){
    chr18.ideo.h2.het <- rbind(chr18.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr19.ideo.h2.het) == 0){
    chr19.ideo.h2.het <- rbind(chr19.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr20.ideo.h2.het) == 0){
    chr20.ideo.h2.het <- rbind(chr20.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr21.ideo.h2.het) == 0){
    chr21.ideo.h2.het <- rbind(chr21.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr22.ideo.h2.het) == 0){
    chr22.ideo.h2.het <- rbind(chr22.ideo.h2.het, no.obs) 
  }
  
  if(nrow(chr1.ideo.h2.hom) == 0){
    chr1.ideo.h2.hom <- rbind(chr1.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr2.ideo.h2.hom) == 0){
    chr2.ideo.h2.hom <- rbind(chr2.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr3.ideo.h2.hom) == 0){
    chr3.ideo.h2.hom <- rbind(chr3.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr4.ideo.h2.hom) == 0){
    chr4.ideo.h2.hom <- rbind(chr4.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr5.ideo.h2.hom) == 0){
    chr5.ideo.h2.hom <- rbind(chr5.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr6.ideo.h2.hom) == 0){
    chr6.ideo.h2.hom <- rbind(chr6.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr7.ideo.h2.hom) == 0){
    chr7.ideo.h2.hom <- rbind(chr7.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr8.ideo.h2.hom) == 0){
    chr8.ideo.h2.hom <- rbind(chr8.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr9.ideo.h2.hom) == 0){
    chr9.ideo.h2.hom <- rbind(chr9.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr10.ideo.h2.hom) == 0){
    chr10.ideo.h2.hom <- rbind(chr10.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr11.ideo.h2.hom) == 0){
    chr11.ideo.h2.hom <- rbind(chr11.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr12.ideo.h2.hom) == 0){
    chr12.ideo.h2.hom <- rbind(chr12.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr13.ideo.h2.hom) == 0){
    chr13.ideo.h2.hom <- rbind(chr13.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr14.ideo.h2.hom) == 0){
    chr14.ideo.h2.hom <- rbind(chr14.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr15.ideo.h2.hom) == 0){
    chr15.ideo.h2.hom <- rbind(chr15.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr16.ideo.h2.hom) == 0){
    chr16.ideo.h2.hom <- rbind(chr16.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr17.ideo.h2.hom) == 0){
    chr17.ideo.h2.hom <- rbind(chr17.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr18.ideo.h2.hom) == 0){
    chr18.ideo.h2.hom <- rbind(chr18.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr19.ideo.h2.hom) == 0){
    chr19.ideo.h2.hom <- rbind(chr19.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr20.ideo.h2.hom) == 0){
    chr20.ideo.h2.hom <- rbind(chr20.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr21.ideo.h2.hom) == 0){
    chr21.ideo.h2.hom <- rbind(chr21.ideo.h2.hom, no.obs) 
  }
  
  if(nrow(chr22.ideo.h2.hom) == 0){
    chr22.ideo.h2.hom <- rbind(chr22.ideo.h2.hom, no.obs) 
  }
  

  if(nrow(chr1.ideo.h2.hom.large) == 0){
    chr1.ideo.h2.hom.large <- rbind(chr1.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr2.ideo.h2.hom.large) == 0){
    chr2.ideo.h2.hom.large <- rbind(chr2.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr3.ideo.h2.hom.large) == 0){
    chr3.ideo.h2.hom.large <- rbind(chr3.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr4.ideo.h2.hom.large) == 0){
    chr4.ideo.h2.hom.large <- rbind(chr4.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr5.ideo.h2.hom.large) == 0){
    chr5.ideo.h2.hom.large <- rbind(chr5.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr6.ideo.h2.hom.large) == 0){
    chr6.ideo.h2.hom.large <- rbind(chr6.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr7.ideo.h2.hom.large) == 0){
    chr7.ideo.h2.hom.large <- rbind(chr7.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr8.ideo.h2.hom.large) == 0){
    chr8.ideo.h2.hom.large <- rbind(chr8.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr9.ideo.h2.hom.large) == 0){
    chr9.ideo.h2.hom.large <- rbind(chr9.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr10.ideo.h2.hom.large) == 0){
    chr10.ideo.h2.hom.large <- rbind(chr10.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr11.ideo.h2.hom.large) == 0){
    chr11.ideo.h2.hom.large <- rbind(chr11.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr12.ideo.h2.hom.large) == 0){
    chr12.ideo.h2.hom.large <- rbind(chr12.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr13.ideo.h2.hom.large) == 0){
    chr13.ideo.h2.hom.large <- rbind(chr13.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr14.ideo.h2.hom.large) == 0){
    chr14.ideo.h2.hom.large <- rbind(chr14.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr15.ideo.h2.hom.large) == 0){
    chr15.ideo.h2.hom.large <- rbind(chr15.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr16.ideo.h2.hom.large) == 0){
    chr16.ideo.h2.hom.large <- rbind(chr16.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr17.ideo.h2.hom.large) == 0){
    chr17.ideo.h2.hom.large <- rbind(chr17.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr18.ideo.h2.hom.large) == 0){
    chr18.ideo.h2.hom.large <- rbind(chr18.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr19.ideo.h2.hom.large) == 0){
    chr19.ideo.h2.hom.large <- rbind(chr19.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr20.ideo.h2.hom.large) == 0){
    chr20.ideo.h2.hom.large <- rbind(chr20.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr21.ideo.h2.hom.large) == 0){
    chr21.ideo.h2.hom.large <- rbind(chr21.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chr22.ideo.h2.hom.large) == 0){
    chr22.ideo.h2.hom.large <- rbind(chr22.ideo.h2.hom.large, no.obs) 
  }
  
  if(nrow(chrx.ideo.haploid.large) == 0){
    chrx.ideo.haploid.large <- rbind(chrx.ideo.haploid.large, no.obs) 
  }
  
  if(nrow(chry.ideo.haploid.large) == 0){
    chry.ideo.haploid.large <- rbind(chry.ideo.haploid.large, no.obs) 
  }
  
  if(nrow(chr1.ideo.h2.het.large) == 0){
    chr1.ideo.h2.het.large <- rbind(chr1.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr2.ideo.h2.het.large) == 0){
    chr2.ideo.h2.het.large <- rbind(chr2.ideo.h2.het.large, no.obs) 
  }
  if(nrow(chr3.ideo.h2.het.large) == 0){
    chr3.ideo.h2.het.large <- rbind(chr3.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr4.ideo.h2.het.large) == 0){
    chr4.ideo.h2.het.large <- rbind(chr4.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr5.ideo.h2.het.large) == 0){
    chr5.ideo.h2.het.large <- rbind(chr5.ideo.h2.het.large, no.obs) 
  }
  if(nrow(chr6.ideo.h2.het.large) == 0){
    chr6.ideo.h2.het.large <- rbind(chr6.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr7.ideo.h2.het.large) == 0){
    chr7.ideo.h2.het.large <- rbind(chr7.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr8.ideo.h2.het.large) == 0){
    chr8.ideo.h2.het.large <- rbind(chr8.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr9.ideo.h2.het.large) == 0){
    chr9.ideo.h2.het.large <- rbind(chr9.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr10.ideo.h2.het.large) == 0){
    chr10.ideo.h2.het.large <- rbind(chr10.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr11.ideo.h2.het.large) == 0){
    chr11.ideo.h2.het.large <- rbind(chr11.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr12.ideo.h2.het.large) == 0){
    chr12.ideo.h2.het.large <- rbind(chr12.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr13.ideo.h2.het.large) == 0){
    chr13.ideo.h2.het.large <- rbind(chr13.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr14.ideo.h2.het.large) == 0){
    chr14.ideo.h2.het.large <- rbind(chr14.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr15.ideo.h2.het.large) == 0){
    chr15.ideo.h2.het.large <- rbind(chr15.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr16.ideo.h2.het.large) == 0){
    chr16.ideo.h2.het.large <- rbind(chr16.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr17.ideo.h2.het.large) == 0){
    chr17.ideo.h2.het.large <- rbind(chr17.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr18.ideo.h2.het.large) == 0){
    chr18.ideo.h2.het.large <- rbind(chr18.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr19.ideo.h2.het.large) == 0){
    chr19.ideo.h2.het.large <- rbind(chr19.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr20.ideo.h2.het.large) == 0){
    chr20.ideo.h2.het.large <- rbind(chr20.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr21.ideo.h2.het.large) == 0){
    chr21.ideo.h2.het.large <- rbind(chr21.ideo.h2.het.large, no.obs) 
  }
  
  if(nrow(chr22.ideo.h2.het.large) == 0){
    chr22.ideo.h2.het.large <- rbind(chr22.ideo.h2.het.large, no.obs) 
  }

  if(nrow(chr1.ideo.h1.het.large) == 0){
    chr1.ideo.h1.het.large <- rbind(chr1.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr2.ideo.h1.het.large) == 0){
    chr2.ideo.h1.het.large <- rbind(chr2.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr3.ideo.h1.het.large) == 0){
    chr3.ideo.h1.het.large <- rbind(chr3.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr4.ideo.h1.het.large) == 0){
    chr4.ideo.h1.het.large <- rbind(chr4.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr5.ideo.h1.het.large) == 0){
    chr5.ideo.h1.het.large <- rbind(chr5.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr6.ideo.h1.het.large) == 0){
    chr6.ideo.h1.het.large <- rbind(chr6.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr7.ideo.h1.het.large) == 0){
    chr7.ideo.h1.het.large <- rbind(chr7.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr8.ideo.h1.het.large) == 0){
    chr8.ideo.h1.het.large <- rbind(chr8.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr9.ideo.h1.het.large) == 0){
    chr9.ideo.h1.het.large <- rbind(chr9.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr10.ideo.h1.het.large) == 0){
    chr10.ideo.h1.het.large <- rbind(chr10.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr11.ideo.h1.het.large) == 0){
    chr11.ideo.h1.het.large <- rbind(chr11.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr12.ideo.h1.het.large) == 0){
    chr12.ideo.h1.het.large <- rbind(chr12.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr13.ideo.h1.het.large) == 0){
    chr13.ideo.h1.het.large <- rbind(chr13.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr14.ideo.h1.het.large) == 0){
    chr14.ideo.h1.het.large <- rbind(chr14.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr15.ideo.h1.het.large) == 0){
    chr15.ideo.h1.het.large <- rbind(chr15.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr16.ideo.h1.het.large) == 0){
    chr16.ideo.h1.het.large <- rbind(chr16.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr17.ideo.h1.het.large) == 0){
    chr17.ideo.h1.het.large <- rbind(chr17.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr18.ideo.h1.het.large) == 0){
    chr18.ideo.h1.het.large <- rbind(chr18.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr19.ideo.h1.het.large) == 0){
    chr19.ideo.h1.het.large <- rbind(chr19.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr20.ideo.h1.het.large) == 0){
    chr20.ideo.h1.het.large <- rbind(chr20.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr21.ideo.h1.het.large) == 0){
    chr21.ideo.h1.het.large <- rbind(chr21.ideo.h1.het.large, no.obs) 
  }
  
  if(nrow(chr22.ideo.h1.het.large) == 0){
    chr22.ideo.h1.het.large <- rbind(chr22.ideo.h1.het.large, no.obs) 
  }

  if(nrow(chr1.ideo.h1.hom.large) == 0){
    chr1.ideo.h1.hom.large <- rbind(chr1.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr2.ideo.h1.hom.large) == 0){
    chr2.ideo.h1.hom.large <- rbind(chr2.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr3.ideo.h1.hom.large) == 0){
    chr3.ideo.h1.hom.large <- rbind(chr3.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr4.ideo.h1.hom.large) == 0){
    chr4.ideo.h1.hom.large <- rbind(chr4.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr5.ideo.h1.hom.large) == 0){
    chr5.ideo.h1.hom.large <- rbind(chr5.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr6.ideo.h1.hom.large) == 0){
    chr6.ideo.h1.hom.large <- rbind(chr6.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr7.ideo.h1.hom.large) == 0){
    chr7.ideo.h1.hom.large <- rbind(chr7.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr8.ideo.h1.hom.large) == 0){
    chr8.ideo.h1.hom.large <- rbind(chr8.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr9.ideo.h1.hom.large) == 0){
    chr9.ideo.h1.hom.large <- rbind(chr9.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr10.ideo.h1.hom.large) == 0){
    chr10.ideo.h1.hom.large <- rbind(chr10.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr11.ideo.h1.hom.large) == 0){
    chr11.ideo.h1.hom.large <- rbind(chr11.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr12.ideo.h1.hom.large) == 0){
    chr12.ideo.h1.hom.large <- rbind(chr12.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr13.ideo.h1.hom.large) == 0){
    chr13.ideo.h1.hom.large <- rbind(chr13.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr14.ideo.h1.hom.large) == 0){
    chr14.ideo.h1.hom.large <- rbind(chr14.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr15.ideo.h1.hom.large) == 0){
    chr15.ideo.h1.hom.large <- rbind(chr15.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr16.ideo.h1.hom.large) == 0){
    chr16.ideo.h1.hom.large <- rbind(chr16.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr17.ideo.h1.hom.large) == 0){
    chr17.ideo.h1.hom.large <- rbind(chr17.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr18.ideo.h1.hom.large) == 0){
    chr18.ideo.h1.hom.large <- rbind(chr18.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr19.ideo.h1.hom.large) == 0){
    chr19.ideo.h1.hom.large <- rbind(chr19.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr20.ideo.h1.hom.large) == 0){
    chr20.ideo.h1.hom.large <- rbind(chr20.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr21.ideo.h1.hom.large) == 0){
    chr21.ideo.h1.hom.large <- rbind(chr21.ideo.h1.hom.large, no.obs) 
  }
  
  if(nrow(chr22.ideo.h1.hom.large) == 0){
    chr22.ideo.h1.hom.large <- rbind(chr22.ideo.h1.hom.large, no.obs) 
  }
  
  ################################## plotting #########################################
  
  chr1.ideo = ggplot() +
    # chr table
    geom_segment(data = chr1.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr1.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr1.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr1.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr1.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr1.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr1.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr1.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr1.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr1.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr1.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr1.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr1.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr1.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr1.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr1.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr1.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr1.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr1.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr1 ", color="black", size = 8)
  
  chr2.ideo = ggplot() +
    # chr table
    geom_segment(data = chr2.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr2.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr2.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr2.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr2.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr2.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr2.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr2.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr2.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr2.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr2.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr2.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr2.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr2.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr2.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr2.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr2.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # centromere
    geom_segment(data = chr2.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr2.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr2.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr2 ", color="black", size = 8)
  
  chr3.ideo = ggplot() +
    # chr table
    geom_segment(data = chr3.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr3.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr3.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr3.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr3.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr3.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr3.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr3.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr3.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr3.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr3.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr3.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr3.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr3.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr3.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr3.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr3.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr3.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr3.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr3 ", color="black", size = 8)
  
  chr4.ideo = ggplot() +
    # chr table
    geom_segment(data = chr4.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr4.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr4.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr4.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr4.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr4.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr4.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr4.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr4.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr4.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr4.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr4.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr4.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr4.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr4.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr4.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr4.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr4.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr4.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr4 ", color="black", size = 8)
  
  chr5.ideo = ggplot() +
    # chr table
    geom_segment(data = chr5.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr5.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr5.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr5.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr5.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr5.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr5.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr5.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr5.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr5.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr5.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr5.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr5.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr5.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr5.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr5.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr5.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr5.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr5.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr5 ", color="black", size = 8)
  
  chr6.ideo = ggplot() +
    # chr table
    geom_segment(data = chr6.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr6.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr6.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr6.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr6.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr6.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr6.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr6.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr6.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr6.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr6.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr6.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr6.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr6.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr6.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr6.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr6.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr6.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr6.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr6 ", color="black", size = 8)
  
  chr7.ideo = ggplot() +
    # chr table
    geom_segment(data = chr7.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr7.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr7.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr7.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr7.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr7.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr7.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr7.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr7.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr7.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr7.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr7.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr7.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr7.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr7.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr7.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr7.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr7.h1.n, color="black", size = 3) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr7.h2.n, color="black", size = 3) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr7 ", color="black", size = 8)
  
  chr8.ideo = ggplot() +
    # chr table
    geom_segment(data = chr8.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr8.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr8.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr8.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr8.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr8.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr8.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr8.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr8.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr8.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr8.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr8.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr8.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr8.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr8.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr8.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr8.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr8.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr8.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr8 ", color="black", size = 8)
  
  chr9.ideo = ggplot() +
    # chr table
    geom_segment(data = chr9.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr9.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr9.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr9.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr9.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr9.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr9.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr9.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr9.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr9.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr9.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr9.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr9.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr9.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr9.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr9.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr9.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr9.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr9.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr9 ", color="black", size = 8)
  
  chr10.ideo = ggplot() +
    # chr table
    geom_segment(data = chr10.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr10.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr10.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr10.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr10.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr10.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr10.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr10.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr10.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr10.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr10.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr10.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr10.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr10.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr10.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr10.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr10.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr10.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr10.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr10", color="black", size = 8)
  
  chr11.ideo = ggplot() +
    # chr table
    geom_segment(data = chr11.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr11.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr11.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr11.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr11.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr11.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr11.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr11.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr11.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr11.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr11.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr11.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr11.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr11.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr11.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr11.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr11.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr11.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr11.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr11", color="black", size = 8)
  
  chr12.ideo = ggplot() +
    # chr table
    geom_segment(data = chr12.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr12.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr12.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr12.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr12.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr12.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr12.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr12.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr12.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr12.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr12.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr12.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr12.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr12.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr12.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr12.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr12.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr12.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr12.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr12", color="black", size = 8)
  
  chr13.ideo = ggplot() +
    # chr table
    geom_segment(data = chr13.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr13.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr13.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr13.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr13.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr13.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr13.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr13.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr13.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr13.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr13.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr13.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr13.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr13.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr13.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr13.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr13.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr13.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr13.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr13", color="black", size = 8)
  
  chr14.ideo = ggplot() +
    # chr table
    geom_segment(data = chr14.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr14.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr14.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr14.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr14.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr14.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr14.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr14.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr14.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr14.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr14.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr14.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr14.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr14.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr14.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr14.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr14.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr14.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr14.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr14", color="black", size = 8)
  
  chr15.ideo = ggplot() +
    # chr table
    geom_segment(data = chr15.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr15.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr15.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr15.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr15.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr15.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr15.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr15.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr15.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr15.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr15.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr15.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr15.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr15.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr15.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr15.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr15.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr15.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr15.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr15", color="black", size = 8)
  
  chr16.ideo = ggplot() +
    # chr table
    geom_segment(data = chr16.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr16.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr16.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr16.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr16.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr16.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr16.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr16.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr16.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr16.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr16.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr16.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr16.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr16.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr16.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr16.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr16.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr16.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr16.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr16", color="black", size = 8)
  
  chr17.ideo = ggplot() +
    # chr table
    geom_segment(data = chr17.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr17.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr17.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr17.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr17.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr17.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr17.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr17.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr17.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr17.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr17.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr17.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr17.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr17.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr17.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr17.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr17.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr17.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr17.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr17", color="black", size = 8)
  
  chr18.ideo = ggplot() +
    # chr table
    geom_segment(data = chr18.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr18.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr18.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr18.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr18.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr18.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr18.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr18.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr18.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr18.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr18.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr18.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr18.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr18.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr18.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr18.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr18.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr18.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr18.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr18", color="black", size = 8)
  
  chr19.ideo = ggplot() +
    # chr table
    geom_segment(data = chr19.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr19.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr19.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr19.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr19.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr19.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr19.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr19.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr19.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr19.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr19.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr19.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr19.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr19.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr19.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr19.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr19.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr19.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr19.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr19", color="black", size = 8)
  
  chr20.ideo = ggplot() +
    # chr table
    geom_segment(data = chr20.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr20.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr20.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr20.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr20.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr20.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr20.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr20.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr20.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr20.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr20.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr20.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr20.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr20.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr20.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr20.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr20.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr20.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr20.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr20", color="black", size = 8)
  
  chr21.ideo = ggplot() +
    # chr table
    geom_segment(data = chr21.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr21.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr21.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr21.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr21.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr21.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr21.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr21.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr21.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr21.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr21.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr21.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr21.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr21.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr21.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr21.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr21.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr21.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr21.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr21", color="black", size = 8)
  
  chr22.ideo = ggplot() +
    # chr table
    geom_segment(data = chr22.table, aes(x = V2, xend = V3, y = 1, yend = 1), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    geom_segment(data = chr22.table, aes(x = V2, xend = V3, y = -1, yend = -1), color = "grey46", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chr22.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 10, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # h1.het
  geom_segment(data = chr22.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr22.ideo.h1.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h1.hom
    geom_segment(data = chr22.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr22.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr22.ideo.h1.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h1.het large
    geom_segment(data = chr22.ideo.h1.het.large, aes(x = start, xend = end, y = 1, yend = 1), color = heterozygous.variants.color, size = 5) +
    # h1.hom large
    geom_segment(data = chr22.ideo.h1.hom.large, aes(x = start, xend = end, y = 1, yend = 1), color = homozygous.varaints.color, size = 5) +
    ######### haplotype 2 #########
  # h2.het
  geom_segment(data = chr22.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chr22.ideo.h2.het, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # h2.hom
    geom_segment(data = chr22.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = 0, yend = genotype), color = homozygous.varaints.color) + 
    geom_segment(data = chr22.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = homozygous.varaints.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    geom_segment(data = chr22.ideo.h2.hom, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = "white", lineend = "round", size = 3, stat = "identity", position = position_dodge()) + 
    # h2.het large
    geom_segment(data = chr22.ideo.h2.het.large, aes(x = start, xend = end, y = -1, yend = -1), color = heterozygous.variants.color, size = 5) +
    # h2.hom large
    geom_segment(data = chr22.ideo.h2.hom.large, aes(x = start, xend = end, y = -1, yend = -1), color = homozygous.varaints.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 1, label=chr22.h1.n, color="black", size = 4) +
    annotate(geom = "text", x = -2000000, y = -1, label=chr22.h2.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="Chr22", color="black", size = 8)
  
  chrx.ideo.plot = ggplot() +
    # chr table
    geom_segment(data = chrx.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chrx.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 5, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # haploid variants
  geom_segment(data = chrx.ideo, aes(x = coordinate, xend = coordinate, y = -1, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chrx.ideo, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # haploid large variants
    geom_segment(data = chrx.ideo.haploid.large, aes(x = start, xend = end, y = 0, yend = 0), color = heterozygous.variants.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 0, label=chrx.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="ChrX ", color="black", size = 8)
  
  chry.ideo.plot = ggplot() +
    # chr table
    geom_segment(data = chry.table, aes(x = V2, xend = V3, y = 0, yend = 0), color = "grey71", lineend = "butt", size = 5, stat = "identity", position = position_dodge()) +
    # centromere
    geom_segment(data = chry.cent, aes(x = V2, xend = V3, y = V1, yend = V1), color = "grey87", size = 5, stat = "identity", position = position_dodge()) +
    ######### haplotype 1 #########
  # haploid variants
  geom_segment(data = chry.ideo, aes(x = coordinate, xend = coordinate, y = -1, yend = genotype), color = heterozygous.variants.color) + 
    geom_segment(data = chry.ideo, aes(x = coordinate, xend = coordinate, y = genotype, yend = genotype), color = heterozygous.variants.color, lineend = "round", size = 4, stat = "identity", position = position_dodge()) + 
    # haploid large variants
    geom_segment(data = chry.ideo.haploid.large, aes(x = start, xend = end, y = 0, yend = 0), color = heterozygous.variants.color, size = 5) +
    # annotate with number of inversion per haplotype
    annotate(geom = "text", x = -2000000, y = 0, label=chry.n, color="black", size = 4) +
    #theme
    theme(axis.title.y = element_text(angle=0, vjust = 0.5, hjust = 1, size = 25), axis.ticks = element_blank(), axis.text.x = element_blank(),  axis.text.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.margin=unit(c(-0.35,0,-0.40,0.05), "null")) + 
    ylab("") + 
    xlab("") +
    ylim(-70, 70) +
    xlim(-15000000, 248956422) +
    annotate(geom = "text", x = -15000000, y = 0, label="ChrY ", color="black", size = 8)
  
    # plot title
  text = paste0(sample.name, "\n GRCh38")
  plot.title = ggplot() + 
    annotate("text", x = 4, y = 25, size=8, label = text) + 
    theme_void()
  
  # export plots
  # ggsave("ideograms.pdf", ideograms, path = "out/", limitsize = FALSE, scale = 1, width = 70, height = 70, units = c("in"), dpi = 300)
  ggsave("plot.title.png", plot.title, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 68.245, height = 3, units = c("cm"), dpi = 300)
  ggsave("chr1.png", chr1.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr2.png", chr2.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr3.png", chr3.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr4.png", chr4.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr5.png", chr5.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr6.png", chr6.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr7.png", chr7.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr8.png", chr8.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr9.png", chr9.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr10.png", chr10.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr11.png", chr11.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr12.png", chr12.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr13.png", chr13.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr14.png", chr14.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr15.png", chr15.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr16.png", chr16.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr17.png", chr17.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr18.png", chr18.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr19.png", chr19.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr20.png", chr20.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr21.png", chr21.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chr22.png", chr22.ideo, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chrx.png", chrx.ideo.plot, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
  ggsave("chry.png", chry.ideo.plot, path = "out/ideograms/", limitsize = FALSE, scale = 1, width = 40, height = 7.8, units = c("cm"), dpi = 300)
