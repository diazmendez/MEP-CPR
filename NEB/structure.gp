unset key
set xrange [-.1:2*pi+1]
set yrange [-.1:2*pi+1]
set isosamples 50
set hidden
set xlabel "\{cita}_1"
set xlabel "\{cita}_2"
set view 30,20,1

set term pos enhanced color

sp -(cos(x)*cos(y))-(sin(x)*sin(y))-(4*((cos(x)**2)+(cos(y)**2))); 
rep 'out' u 1:2:((-(cos($1)*cos($2))-(sin($1)*sin($2))-(4*((cos($1)**2)+(cos($2)**2))))+1)  w lp

set output "surface.eps"
rep

