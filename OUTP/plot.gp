set xlabel "nump"
set ylabel "time"
set terminal pdf
set output "times.pdf"
plot \
"../times.csv" using 1:3 with linespoints title "compF", \
"" u 1:4 w lp  title "compV"
