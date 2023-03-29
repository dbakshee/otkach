ne = 80*3
nf = 200*3
# 1001
set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size nf, ne
unset key
unset colorbox
set view map

bin = "w02/w02_fftdos_pha.bin"
png = "w02/w02_fftdos_pha.png"
set output png
plot bin matrix binary with image

bin = "w02/w02_fftdos_amp.bin"
png = "w02/w02_fftdos_amp.png"
set output png
plot bin matrix binary with image
