set term postscript dashed dl 4 lw 1.2 eps enhanced color 27 rounded
set size 0.9,1.0

# Useful physical and mathematical formulas:

hbc = 197.33
Pi = 3.14159265

set style function filledcurves above y1=0.0
set style fill transparent solid 0.5 noborder

set datafile separator ','

###########################
###########################
#### BW distributions: #### 
###########################
###########################

set output 'M12-fit12-BW-bins.eps'

set xlabel "M(J/{/Symbol y} p)" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [4.03:4.43][0:70] 'M12-fit12-BW-histo.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory (BW)", 'Fig2e.csv' u 4:5:6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):5:6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}"

set output 'Mz-fit12-BW-bins.eps'

set xlabel "cos({/Symbol q}_p)" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-BW-histo.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory (BW)" , 'Fig2c.csv' u 4:6:7 w e lw 3 lc 2 notitle

set output 'M12-fit23-BW-bins.eps'

set xlabel "M(J/{/Symbol y} p)" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [4.03:4.43][0:70] 'M12-fit23-BW-histo.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory (BW)", 'Fig2e.csv' u 4:5:6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):5:6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}"

set output 'Mz-fit23-BW-bins.eps'

set xlabel "cos({/Symbol q}_p)" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit23-BW-histo.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory (BW)" , 'Fig2c.csv' u 4:6:7 w e lw 3 lc 2 notitle

#####################################
#####################################
#### Rescattering distributions: ####
#####################################
#####################################

################
# Scenario (a) #
################

set output 'M12-fit12-T4-bins-a.eps'

set xlabel "M(J/{/Symbol y} p) [a]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-histo-a.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2e.csv' u 4:5:6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):5:6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}"

set output 'Mz-fit12-T4-bins-a.eps'

set xlabel "cos({/Symbol q}_p) [a]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-histo-a.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2c.csv' u 4:6:7 w e lw 3 lc 2 t "J/{/Symbol y} p"

set output 'M12-fit12-T4-a.eps'

set xlabel "M(J/{/Symbol y} p) [a]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-a.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2e.csv' u 4:($5-$7):6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):($5-$7):6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}", 'M12-fit12-T4-PD-a.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle


set output 'Mz-fit12-T4-a.eps'

set xlabel "cos({/Symbol q}_p) [a]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-a.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2c.csv' u 4:($6-$8):7 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Mz-fit12-T4-PD-a.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle

################
# Scenario (b) #
################

set output 'M12-fit12-T4-bins-b.eps'

set xlabel "M(J/{/Symbol y} p) [b]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-histo-b.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2e.csv' u 4:5:6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):5:6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}"

set output 'Mz-fit12-T4-bins-b.eps'

set xlabel "cos({/Symbol q}_p) [b]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-histo-b.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2c.csv' u 4:6:7 w e lw 3 lc 2 t "J/{/Symbol y} p"

set output 'M12-fit12-T4-b.eps'

set xlabel "M(J/{/Symbol y} p) [b]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-b.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2e.csv' u 4:($5-$7):6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):($5-$7):6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}", 'M12-fit12-T4-PD-b.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle


set output 'Mz-fit12-T4-b.eps'

set xlabel "cos({/Symbol q}_p) [b]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-b.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2c.csv' u 4:($6-$8):7 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Mz-fit12-T4-PD-b.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle


################
# Scenario (c) #
################

set output 'M12-fit12-T4-bins-c.eps'

set xlabel "M(J/{/Symbol y} p) [c]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-histo-c.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2e.csv' u 4:5:6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):5:6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}"

set output 'Mz-fit12-T4-bins-c.eps'

set xlabel "cos({/Symbol q}_p) [c]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-histo-c.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2c.csv' u 4:6:7 w e lw 3 lc 2 t "J/{/Symbol y} p"

set output 'M12-fit12-T4-c.eps'

set xlabel "M(J/{/Symbol y} p) [c]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-c.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2e.csv' u 4:($5-$7):6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):($5-$7):6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}", 'M12-fit12-T4-PD-c.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle


set output 'Mz-fit12-T4-c.eps'

set xlabel "cos({/Symbol q}_p) [c]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-c.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2c.csv' u 4:($6-$8):7 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Mz-fit12-T4-PD-c.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle


################
# Scenario (d) #
################

set output 'M12-fit12-T4-bins-d.eps'

set xlabel "M(J/{/Symbol y} p) [d]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-histo-d.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2e.csv' u 4:5:6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):5:6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}"

set output 'Mz-fit12-T4-bins-d.eps'

set xlabel "cos({/Symbol q}_p) [d]" offset 0,-.5 font "Helvetica,30"
set ylabel "events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-histo-d.plot' u 1:2 w p ps 3 pt 5 lc 6 t "theory", 'Fig2c.csv' u 4:6:7 w e lw 3 lc 2 t "J/{/Symbol y} p"

set output 'M12-fit12-T4-d.eps'

set xlabel "M(J/{/Symbol y} p) [d]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.1
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

set arrow from 4.1513, 0.0 to 4.1513, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.2933, 0.0 to 4.2933, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3175, 0.0 to 4.3175, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)
set arrow from 4.3236, 0.0 to 4.3236, 70 nohead lw 7 lc 7 dt (0.25,0.5,0.25,0.5)

plot [4.03:4.43][0:70] 'M12-fit12-T4-PD-d.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2e.csv' u 4:($5-$7):6 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Fig2f.csv' u ($4+0.002):($5-$7):6 w e lw 2 lc 4 t "J/{/Symbol y} ~p{.4\261}", 'M12-fit12-T4-PD-d.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle


set output 'Mz-fit12-T4-d.eps'

set xlabel "cos({/Symbol q}_p) [d]" offset 0,-.5 font "Helvetica,30"
set ylabel "lineshapes + events / bin" offset -0.5,-0.3 font "Helvetica,30"

set xtics 0.5
set ytics 10.0

set key samplen 0
set key width -4
set key left Left reverse

set errorbars 2

plot [-1:1][0:70] 'Mz-fit12-T4-PD-d.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) t "theory", 'Fig2c.csv' u 4:($6-$8):7 w e lw 3 lc 2 t "J/{/Symbol y} p", 'Mz-fit12-T4-PD-d.plot' u 1:2 w l lw 6 lc 6 dt (1.25,0.5,1.25,0.5) notitle

