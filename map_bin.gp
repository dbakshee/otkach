# W03: 1201x125 - 1920x800
# W05: 1001x181 - 1600x1158
# W02: 2001x47 - 1600x150
# W02: 2001x80 - 1600x255
# W02: 2001x143 - 1600x456

vals = "dos rxy21 rxx rxy12 r2t22 r2t11 ryy"

set terminal pngcairo  notransparent enhanced font "arial,10" fontscale 1.0 size 1920, 800
unset key
# unset colorbox
set xrange noextend
set xtics out
set yrange noextend
set ytics out

set view map
do for [val in "dos"] {
    bin = "w02/w02a_".val.".bin"
    png = "w02/w02a_".val.".png"
    set output png
    plot bin matrix binary with image
}
