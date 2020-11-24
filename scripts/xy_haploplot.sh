#!/bin/bash

Rscript ./scripts/xy_haploplot.r
Rscript ./scripts/xy_bed_reads_sub.R
sh ./scripts/xy_im.crop.sh
sh ./scripts/bedtools.sub.sh
perl ./scripts/pdfplot/scripts/urls2pdf_hgsid.pl ./out/bed/*.bed ./out/ideograms/ideogram.pdf male
rm out/ideograms/ideogram.pdf
