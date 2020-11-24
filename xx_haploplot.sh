#!/bin/bash

Rscript ./scripts/xx_haploplot.r
sh ./scripts/xx_im.crop.sh
hgsid=$(perl ./scripts/pdfplot/scripts/urls2pdf_hgsid.pl ./out/bed/*.bed ./out/ideograms/ideogram.pdf female)
rm out/ideograms/ideogram.pdf


for bed in $(ls ./in/reads_bed/*.bed ./in/reads_bed/*.bed.gz)
do

perl ./scripts/pdfplot/scripts/upload.pl "$bed" "$hgsid"

done
