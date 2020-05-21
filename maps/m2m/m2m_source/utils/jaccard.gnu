set term post eps enh color
set out 'jaccard.eps'
set size ratio -1

set xlabel 'Cell x-index'
set ylabel 'Cell y-index'
#==============================================

#set palette maxcolors 5
set palette defined ( 0 '#000fff',\
                      1 '#90ff70',\
                      2 '#ee0000')

#set palette rgbformulae 33,13,10

# To set the scale of colored values.
set cbrange [0:1]
set zrange [0:1]

#set view map

#NOTE: to use this option there must be blank lines in data file every time x-coordinate changes.
set pm3d map

# We use awk utility to reformat the data.

splot '< awk -f addblanks.awk ./output/jaccard_indexes.txt' u 1:2:8 w pm3d palette notitle

