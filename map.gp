set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 1920, 800
set output 'dos.png'
unset key
unset colorbox
set view map
plot 'dos.dat' matrix with image
