#!/bin/bash

Rscript ./scripts/xy_haploplot.r
Rscript ./scripts/bed_reads_sub.r
sh ./scripts/xy_im.crop.sh
sh ./scripts/bedtools.sub.sh
hgsid=$(perl ./scripts/pdfplot/scripts/urls2pdf_hgsid.pl ./out/bed/*.bed ./out/ideograms/ideogram.pdf male)
rm out/ideograms/ideogram.pdf


for bed in $(ls ./out/bed_reads/w?.c?.adj.bed.gz)
do

perl ./scripts/pdfplot/scripts/upload.pl "$bed" "$hgsid"

done
