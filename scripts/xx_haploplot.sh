#!/bin/bash

Rscript ./scripts/xx_haploplot.r
sh ./scripts/xx_im.crop.sh
perl ./scripts/pdfplot/scripts/urls2pdf_hgsid.pl ./out/bed/*.bed ./out/ideograms/ideogram.pdf female
rm out/ideograms/ideogram.pdf
