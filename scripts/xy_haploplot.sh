#!/bin/bash

Rscript ./scripts/xy_haploplot.r
sh ./scripts/xy_im.crop.sh
Rscript ./scripts/bed_reads_sub.r
sh ./scripts/bedtools.sub.sh
$hgsid=(perl ./scripts/pdfplot/scripts/urls2pdf_hgsid.pl ./out/bed/*.bed ./out/ideograms/ideogram.pdf male)
rm out/ideograms/ideogram.pdf



for bed in $(ls ./in/reads_bed/*.bed ./in/reads_bed/*.bed.gz)
do

perl ./scripts/pdfplot/scripts/upload.pl "$bed" "$hgsid"

done

