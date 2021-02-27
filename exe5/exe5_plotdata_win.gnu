#!/usr/bin/gnuplot
#
# Plotting a color map using the default Matlab palette
#
# AUTHOR: Gill
# REFERENCE: 
#	Hagen Wierstorf (Original Script)
#	http://www.gnuplotting.org/code/default_color_map3.gnu
#   http://www.gnuplotting.org/tag/colormap/

reset

a = ARG1
b = ARG2
dx = ARG3
dy = ARG4
u = ARG5
v = ARG6
tmax = ARG7
file = ARG8
dir = ARG9

# wxt
set terminal wxt size 1050,720 enhanced font 'Verdana,10' persist
# png
set terminal pngcairo size 1024,720 enhanced font 'Verdana,12'
print dir
set output dir."\\exe5_plot_colormap.png"

unset key

# border
set style line 11 lc rgb '#808080' lt 1
set border 0 front ls 11
set tics nomirror out scale 1
set xrange [0:a]
set yrange [0:b]
set xtics 0,dx,a
set ytics 0,dy,b
set title "SOLUÇÃO NUMÉRICA POR VOLUMES FINITOS\n [dx=".dx.", dy=".dy.", u=".u.", v=".v."]"
set xlabel 'a(m)'	
set ylabel 'b(m)'
# disable colorbar tics
set cbtics scale 0
set cbrange [0:1]
set cblabel "Temperatura (K)"
# matlab palette colors
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")

plot file u ($1+dx/2):($2+dy/2):3 w image

#INTERPOLAÇÃO
set output dir."_exe5_plot_color_int.png"
set pm3d map
set pm3d interpolate 0,0
splot file u ($1+dx/2):($2+dy/2):3

#LINHAS
set output dir."_exe5_plot_lines.png"
#plot file u ($1+dx/2):($2+dy/2):($3) w linespoints ls 1
set ylabel 'T(k)'
set xrange [0:a]
set ytics 0,tmax*0.1,(tmax+tmax*0.1)
set yrange [0:(tmax+tmax*0.1)]
plot file u ($1+dx/2):3 w linespoints ls 1