set terminal pngcairo  enhanced font "arial,10" fontscale 1.0 size 1920, 600
o=0.002
if (1) {
    set output 'dos_t.png'
    plot for [i=2:50] 'dos_t.dat' using 1:(column(i) + (i-2)*o) title sprintf("%.1f", 8.0+(i-2.0)/10) with lines
} else {
    set output 'dos_tb.png'
    plot for [i=2:50] 'dos_t.dat' using (1/$1):(column(i) + (i-2)*o) title sprintf("%.1f", 8.0+(i-2.0)/10) with lines
}

