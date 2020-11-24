#!/bin/bash
bedtools intersect -a in/bed_reads/*.WW.CC.bam_reads.bed.gz -b in/bed_reads/2x.adjascent.bed -u > out/bed_reads/ww.cc.adj.bed
sed -i '1s/^/track name=WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' out/bed_reads/ww.cc.adj.bed

bedtools intersect -a in/bed_reads/*.WC.CW.bam_reads.bed.gz -b in/bed_reads/2x.adjascent.bed -u > out/bed_reads/wc.cw.adj.bed
sed -i '1s/^/track name=WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' out/bed_reads/wc.cw.adj.bed

gzip out/bed_reads/*.bed
