#!/bin/bash

Rscript ./scripts/xy_haploplot.r
sh ./scripts/xy_im.crop.sh
perl ./scripts/pdfplot/scripts/urls2pdf_hgsid.pl ./out/bed/*.bed ./out/ideograms/ideogram.pdf male
rm out/ideograms/ideogram.pdf
