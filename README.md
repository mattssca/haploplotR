  # haploplotR
R package for visualizing haplotype-resolved inversions, associated metrics and variant characteristics

* This package takes output produced by InvertypeR and splits data table into two data frames depending on genotype (H1 and H2)
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
5. Execute `haploplot_run.sh`

* To add hyperlinks from UCSC genome browser to ideograms:
  * Execute bedman.R scripts (scripts/pdfplot/scripts/)
  * Execute urls2pdf_hgsid.pl script (scripts/pdfplot/scripts/) to parse variants with UCSC Genome Browser links

### Dependencies
Package | Version | Enviroment
--------| --------|-----------
dplyr | 0.8.5 | R
gridExtra | 2.3 | R
ggplot2 | 3.3.0 | R
data.table |1.12.8 | R
psych |2.0.8 | R
ImageMagick | 7.0.10-31 | bash
