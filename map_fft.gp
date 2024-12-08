ne = 80*3
nf = 200*3
# 1001
set terminal pngcairo  notransparent enhanced font "arial,10" fontscale 1.0 size 600, 240
unset key
#unset colorbox
set view map

set xrange [0:*] noextend
set xtics out
set yrange noextend
set ytics out
set ylabel "Energy"
set xlabel "Frequency (per unit of ob)"

wxx = "w02"
bin = wxx."/".wxx."a_fftdos_pha.bin"
png = wxx."/".wxx."a_fftdos_pha.png"
set output png
set title "Phase of RFFT(DOS) w=0.2"
plot bin matrix binary with image

bin = wxx."/".wxx."a_fftdos_amp.bin"
png = wxx."/".wxx."a_fftdos_amp.png"
set output png
set title "Log Amp of RFFT(DOS) w=0.2"
plot bin matrix binary with image
