  # haploplotR
Program for visualizing haplotype-resolved inversions, associated metrics and variant characteristics

* This set of scripts takes output produced by InvertypeR and splits data table into two data frames depending on genotype (H1 and H2)
* Haplotypes are filtered on blacklisted regions and variants intersecting centromeres
* Package outputs:
  * Clickable haplotype-resolved ideograms annotating heterozygous and homozygous variants along each chromosome
  * Genome-wide bed-files (H1 and H2) formatted for genome browser interrogation
  * Figures related to variant-size distribution
  * Variant summary metrics table
  * Haplotype speciic tables with annotations such as; probabillity scores, read counts (WW, WC, CC)
  
### Intructions
1. Install dependencies
2. Clone/download repository
2. cd main folder
4. Put output from InvertypeR in Input folder
5. Put composite read-data browser files (BPR output) in  in/bed_reads/ (wc.cw.bed.gz and ww.cc.bed.gz)
5. If session ID for UCSC genome browser is avaialble, add ID at end of line 5 in scripts/nn_haploplot.sh script, if not unique session ID will be generated
6. Execute `haploplot_run.sh` and follow prompted instructions

### Dependencies
Package | Version | Enviroment
--------| --------|-----------
dplyr | 0.8.5 | R
gridExtra | 2.3 | R
ggplot2 | 3.3.0 | R
data.table |1.12.8 | R
psych |2.0.8 | R
bedtools |2.26 | bash
[ImageMagick](https://imagemagick.org/index.php) | 7.0.10-31 | bash
[img2pdf](https://pypi.org/project/img2pdf/) | 0.4.0 | python
[PDF-API2](https://metacpan.org/pod/PDF::API2::Simple) | 1.1.14.u | perl
[LBW::UserAgent](https://metacpan.org/pod/LWP::UserAgent) | 4.69 | perl

