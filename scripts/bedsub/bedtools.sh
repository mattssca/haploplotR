#!/bin/bash

##########################################################HG005#####################################################################
bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00512.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.hg00512.bed -u > Desktop/NEW\ BED/READS/filtered/hg00512.ww.cc.adj.bed
sed -i '1s/^/track name=HG00512.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/hg00512.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/hg00512.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00512.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.hg00512.bed -u > Desktop/NEW\ BED/READS/filtered/hg00512.wc.cw.adj.bed
sed -i '1s/^/track name=HG00512.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/hg00512.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/hg00512.wc.cw.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00513.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.hg00513.bed -u > Desktop/NEW\ BED/READS/filtered/hg00513.ww.cc.adj.bed
sed -i '1s/^/track name=HG00513.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/hg00513.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/hg00513.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00513.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.hg00513.bed -u > Desktop/NEW\ BED/READS/filtered/hg00513.wc.cw.adj.bed
sed -i '1s/^/track name=HG00513.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/hg00514.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/hg00513.wc.cw.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00514.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.hg00514.bed -u > Desktop/NEW\ BED/READS/filtered/hg00514.ww.cc.adj.bed
sed -i '1s/^/track name=HG00514.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/hg00514.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/hg00514.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00514.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.hg00514.bed -u > Desktop/NEW\ BED/READS/filtered/hg00514.wc.cw.adj.bed
sed -i '1s/^/track name=HG00514.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/hg00514.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/hg00514.wc.cw.adj.bed

sed -i '1s/^/track name=HG00512.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/hg00512.txt
sed -i '1s/^/track name=HG00513.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/hg00513.txt
sed -i '1s/^/track name=HG00514.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/hg00514.txt

##########################################################HG007#####################################################################
bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00731.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.HG00731.bed -u > Desktop/NEW\ BED/READS/filtered/HG00731.ww.cc.adj.bed
sed -i '1s/^/track name=HG00731.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/HG00731.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/HG00731.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00731.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.HG00731.bed -u > Desktop/NEW\ BED/READS/filtered/HG00731.wc.cw.adj.bed
sed -i '1s/^/track name=HG00731.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/HG00731.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/HG00731.wc.cw.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00732.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.HG00732.bed -u > Desktop/NEW\ BED/READS/filtered/HG00732.ww.cc.adj.bed
sed -i '1s/^/track name=HG00732.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/HG00732.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/HG00732.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00732.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.HG00732.bed -u > Desktop/NEW\ BED/READS/filtered/HG00732.wc.cw.adj.bed
sed -i '1s/^/track name=HG00732.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/HG00733.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/HG00732.wc.cw.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00733.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.HG00733.bed -u > Desktop/NEW\ BED/READS/filtered/HG00733.ww.cc.adj.bed
sed -i '1s/^/track name=HG00733.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/HG00733.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/HG00733.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/HG00733.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.HG00733.bed -u > Desktop/NEW\ BED/READS/filtered/HG00733.wc.cw.adj.bed
sed -i '1s/^/track name=HG00733.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/HG00733.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/HG00733.wc.cw.adj.bed

sed -i '1s/^/track name=HG00731.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/hg00731.txt
sed -i '1s/^/track name=HG00732.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/hg00732.txt
sed -i '1s/^/track name=HG00733.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/hg00733.txt

##########################################################NA192#####################################################################
bedtools intersect -a Desktop/NEW\ BED/READS/full/NA19238.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.NA19238.bed -u > Desktop/NEW\ BED/READS/filtered/NA19238.ww.cc.adj.bed
sed -i '1s/^/track name=NA19238.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/NA19238.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/NA19238.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/NA19238.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.NA19238.bed -u > Desktop/NEW\ BED/READS/filtered/NA19238.wc.cw.adj.bed
sed -i '1s/^/track name=NA19238.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/NA19238.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/NA19238.wc.cw.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/NA19239.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.NA19239.bed -u > Desktop/NEW\ BED/READS/filtered/NA19239.ww.cc.adj.bed
sed -i '1s/^/track name=NA19239.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/NA19239.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/NA19239.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/NA19239.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.NA19239.bed -u > Desktop/NEW\ BED/READS/filtered/NA19239.wc.cw.adj.bed
sed -i '1s/^/track name=NA19239.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/NA19240.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/NA19239.wc.cw.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/NA19240.WW.CC.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.NA19240.bed -u > Desktop/NEW\ BED/READS/filtered/NA19240.ww.cc.adj.bed
sed -i '1s/^/track name=NA19240.WW.CC.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/NA19240.ww.cc.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/NA19240.ww.cc.adj.bed

bedtools intersect -a Desktop/NEW\ BED/READS/full/NA19240.WC.CW.bam_reads.bed -b Desktop/NEW\ BED/adjascent.regions/adj.NA19240.bed -u > Desktop/NEW\ BED/READS/filtered/NA19240.wc.cw.adj.bed
sed -i '1s/^/track name=NA19240.WC.CW.reads visibility=1 colorByStrand="103,139,139 243,165,97"\n/' Desktop/NEW\ BED/READS/filtered/NA19240.wc.cw.adj.bed
gzip Desktop/NEW\ BED/READS/filtered/NA19240.wc.cw.adj.bed

sed -i '1s/^/track name=NA19238.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/na19238.txt
sed -i '1s/^/track name=NA19239.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/na19239.txt
sed -i '1s/^/track name=NA19240.inversions color="52, 115, 37"\n/' Desktop/NEW\ BED/VARIANTS/na19240.txt