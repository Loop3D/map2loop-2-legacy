# This script inserts a blank line every time the value in specified column changes. (Need for pm3d in gnuplot).
(NR <= 1) {next} # skip first line (header)
(NR == 2) {prev = $3}
($3 != prev) {print ""; prev = $3} # print blank line  
{print} # print the line
