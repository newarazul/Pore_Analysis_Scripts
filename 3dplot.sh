set terminal postscript eps enhanced color size 4in,3in
set output 'TRy.eps'

set style data lines
set contour base
set surface
set pm3d
set xrange [-7:7]
set yrange [-7:7]
set zrange [0:2]
set view 30,10,1.0,1.0
set dgrid3d 40,51,3




splot "flat.dat"
