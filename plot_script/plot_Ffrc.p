# Gnuplot script for plotting data
set terminal epslatex color
set out '~/PhD/code/C_Finite_Temp/pics/Ffrc.tex'
set arrow from 1E-3,1E-14 to 1E-3,1E-14 lw 5 lc rgb 'black'
set xrange [1.E-4 : 1E-2]
set yrange [5.E-13 : 1.E-1]
set key Left center top reverse
set key at 1.9E-3,1E-2
set key font ",8"
set key spacing 4.5
set logscale x 10
set logscale y 10
set format x '$10^{%T}$'
set format y '$10^{%T}$'
set border lw 3 back
set size 1,0.5
set ytics 1e-12, 1e4, 1e-2 offset 0.5,0
set xtics 1e-4, 1e1, 1e-2
set xlabel "$v/c$" offset 0,.5
set ylabel '$F~[\mathrm{fN}]$' offset 3,0
set label at graph 0.8E0,0.3 '\Large$\propto\frac{v^3}{(2z_a)^{10}}$'

fil = "../output/resultsfix.dat"
hbar = 1.0545718E-34
wsp1 = 0.8/6.5821220E-16
a0   = 5.26E-39
eps0 = 8.854187817E-12
c    = 299792458.
pi   = 3.1415926535897932384626433
    
F0 = 3*hbar*wsp1**5*a0/(2*pi*eps0*c**4)*1E15


set multiplot

plot fil u ($1):($2*F0) smooth unique lw 4 dashtype 2 lc rgb 'red' title '$F^\mathrm{t}$',\
     fil u ($1):(-$3*F0) smooth unique lw 4 dashtype 3 lc rgb 'black' title '$F^\mathrm{r}$',\
     fil u ($1):(($2+$3)*F0) smooth unique lw 4 lc rgb 'blue' notitle,\
     fil u ($1):(($4+$5)*F0) smooth unique lw 3 dashtype 4 lc rgb 'black' notitle

set origin 0.14,0.17
set size 0.6,0.3
unset ylabel
unset xlabel
unset logscale y
unset label
set format y2 '\small %g'
set format x '\small $10^{%T}$'
set yrange [0. : 1.]
set y2tics 0, 0.3, 1 offset -0.5
set xtics 1e-4, 1e1, 0.9e-2 offset 1,0.4
unset key

set label at graph 0.4E0,0.45 '$\left|F^\mathrm{r}/F^\mathrm{t}\right|$'
unset ytics


plot fil u ($1):(-$3/$2) smooth unique lw 4 lc rgb 'black' notitle

unset multiplot

     
set out
