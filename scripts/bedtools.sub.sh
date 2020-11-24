#!/bin/bash
bedtools intersect -a in/bed_reads/*.WW.CC.bam_reads.bed.gz -b 2x.adjascent.bed -u > ww.cc.adj.bed
sed -i '1s/^/track name=WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' h

bedtools intersect -a in/bed_reads/*.WC.CW.bam_reads.bed.gz -b 2x.adjascent.bed -u > wc.cw.adj.bed
sed -i '1s/^/track name=HG00512.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/'

gzip in/bed_reads/*.bed