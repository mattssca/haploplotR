#!/bin/bash

# crop images
convert out/ideograms/chr1.png -crop 4424x193+132+415 out/ideograms/chr1.png
convert out/ideograms/chr2.png -crop 4424x193+132+415 out/ideograms/chr2.png
convert out/ideograms/chr3.png -crop 4424x193+132+415 out/ideograms/chr3.png
convert out/ideograms/chr4.png -crop 4424x193+132+415 out/ideograms/chr4.png
convert out/ideograms/chr5.png -crop 4424x193+132+415 out/ideograms/chr5.png
convert out/ideograms/chr6.png -crop 4424x193+132+415 out/ideograms/chr6.png
convert out/ideograms/chr7.png -crop 4424x193+132+415 out/ideograms/chr7.png
convert out/ideograms/chr8.png -crop 4424x193+132+415 out/ideograms/chr8.png
convert out/ideograms/chr9.png -crop 4424x193+132+415 out/ideograms/chr9.png
convert out/ideograms/chr10.png -crop 4424x193+132+415 out/ideograms/chr10.png
convert out/ideograms/chr11.png -crop 4424x193+132+415 out/ideograms/chr11.png
convert out/ideograms/chr12.png -crop 4424x193+132+415 out/ideograms/chr12.png
convert out/ideograms/chr13.png -crop 3036x193+132+415 out/ideograms/chr13.png
convert out/ideograms/chr14.png -crop 3036x193+132+415 out/ideograms/chr14.png
convert out/ideograms/chr15.png -crop 3036x193+132+415 out/ideograms/chr15.png
convert out/ideograms/chr16.png -crop 3036x193+132+415 out/ideograms/chr16.png
convert out/ideograms/chr17.png -crop 3036x193+132+415 out/ideograms/chr17.png
convert out/ideograms/chr18.png -crop 3036x193+132+415 out/ideograms/chr18.png
convert out/ideograms/chr19.png -crop 3036x193+132+415 out/ideograms/chr19.png
convert out/ideograms/chr20.png -crop 3036x193+132+415 out/ideograms/chr20.png
convert out/ideograms/chr21.png -crop 3036x193+132+415 out/ideograms/chr21.png
convert out/ideograms/chr22.png -crop 3036x193+132+415 out/ideograms/chr22.png
convert out/ideograms/chrx.png -crop 3036x193+132+415 out/ideograms/chrx.png
convert out/ideograms/chry.png -crop 3036x193+132+415 out/ideograms/chry.png
convert out/ideograms/plot.title.png -crop 7460x354+300+0 out/ideograms/plot.title.png 


# compile ideograms fro chr 1-12
convert out/ideograms/chr1.png out/ideograms/chr2.png out/ideograms/chr3.png out/ideograms/chr4.png out/ideograms/chr5.png out/ideograms/chr6.png out/ideograms/chr7.png out/ideograms/chr8.png out/ideograms/chr9.png out/ideograms/chr10.png out/ideograms/chr11.png out/ideograms/chr12.png -append out/ideograms/1-12.ideogram.png

# compile ideograms for chr 13-x(y)
convert out/ideograms/chr13.png out/ideograms/chr14.png out/ideograms/chr15.png out/ideograms/chr16.png out/ideograms/chr17.png out/ideograms/chr18.png out/ideograms/chr19.png out/ideograms/chr20.png out/ideograms/chr21.png out/ideograms/chr22.png out/ideograms/chrx.png out/ideograms/chry.png -append out/ideograms/13-xy.ideogram.png

# add both ideograms horizonatally
convert +append out/ideograms/1-12.ideogram.png out/ideograms/13-xy.ideogram.png out/ideograms/ideogram.png

# add plot title
convert out/ideograms/plot.title.png out/ideograms/ideogram.png	-append out/ideograms/ideogram.png

# clean enviroment and migrate files
convert +append out/figs/binned.box.png out/figs/viol.png out/figs/variantsize.png
rm out/figs/binned.box.png
rm out/figs/viol.png
rm out/ideograms/plot.title.png
rm out/ideograms/1-12.ideogram.png
rm out/ideograms/13-xy.ideogram.png
mv out/ideograms/chr1.png out/ideograms/chr/chr1.png 
mv out/ideograms/chr2.png out/ideograms/chr/chr2.png 
mv out/ideograms/chr3.png out/ideograms/chr/chr3.png 
mv out/ideograms/chr4.png out/ideograms/chr/chr4.png 
mv out/ideograms/chr5.png out/ideograms/chr/chr5.png 
mv out/ideograms/chr6.png out/ideograms/chr/chr6.png 
mv out/ideograms/chr7.png out/ideograms/chr/chr7.png 
mv out/ideograms/chr8.png out/ideograms/chr/chr8.png 
mv out/ideograms/chr9.png out/ideograms/chr/chr9.png 
mv out/ideograms/chr10.png out/ideograms/chr/chr10.png
mv out/ideograms/chr11.png out/ideograms/chr/chr11.png
mv out/ideograms/chr12.png out/ideograms/chr/chr12.png
mv out/ideograms/chr13.png out/ideograms/chr/chr13.png
mv out/ideograms/chr14.png out/ideograms/chr/chr14.png
mv out/ideograms/chr15.png out/ideograms/chr/chr15.png
mv out/ideograms/chr16.png out/ideograms/chr/chr16.png
mv out/ideograms/chr17.png out/ideograms/chr/chr17.png
mv out/ideograms/chr18.png out/ideograms/chr/chr18.png
mv out/ideograms/chr19.png out/ideograms/chr/chr19.png
mv out/ideograms/chr20.png out/ideograms/chr/chr20.png
mv out/ideograms/chr21.png out/ideograms/chr/chr21.png
mv out/ideograms/chr22.png out/ideograms/chr/chr22.png
mv out/ideograms/chrx.png out/ideograms/chr/chrx.png 
mv out/ideograms/chry.png out/ideograms/chr/chry.png 

img2pdf out/ideograms/ideogram.png -o out/ideograms/ideogram.pdf
