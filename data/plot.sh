#!/bin/bash

# preparation: delete first line
# for el in data.csv.*;do sed -i -e "1d" $el; done

# plotting
# mkdir png
for i in `seq 1 1001`;do
gnuplot << EOF
set datafile separator ","
unset key
set terminal png size 800,600
set out "png/frame_${i}.png"
plot "data.csv.${i}" using 3:6 with points
# plot "data.csv.${i}" using 1:2:3:4 with vectors
EOF
done

# video encoding
# mencoder mf://frame_%d.png -ovc copy -o video.avi
